#include "core.hpp"
#include <QDir>
#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QStandardPaths>
#include <QTextStream>
#include <cmath>
#include <functional>
#include <orf/orf.hpp>
namespace go {
static QStringList split(QString line, QChar d) {
    QStringList r;
    QString cell;
    bool quote = false;
    for (int i = 0; i < line.size(); ++i) {
        auto c = line[i];
        if (c == '"') {
            if (quote && i + 1 < line.size() && line[i + 1] == '"') {
                cell += '"';
                ++i;
            } else
                quote = !quote;
        } else if (c == d && !quote) {
            r << cell;
            cell.clear();
        } else
            cell += c;
    }
    if (quote)
        throw homology::Error("Незакрытые кавычки в таблице");
    r << cell;
    return r;
}
QVector<Row> table(const QString &path, QChar d) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        throw homology::Error("Не удалось открыть таблицу: " + path);
    QTextStream in(&f);
    auto names = split(in.readLine().trimmed(), d);
    if (!names.isEmpty())
        names[0].remove(QChar(0xfeff));
    QVector<Row> out;
    int line = 1;
    while (!in.atEnd()) {
        auto l = in.readLine();
        ++line;
        if (l.trimmed().isEmpty())
            continue;
        auto values = split(l, d);
        if (names.size() != values.size())
            throw homology::Error(
                QString("Строка %1 таблицы %2: неверное число столбцов").arg(line).arg(path));
        Row row;
        for (int i = 0; i < names.size(); ++i)
            row[names[i]] = values[i];
        out << row;
    }
    return out;
}
void writeTable(const QString &path, const QStringList &columns, const QVector<Row> &rows) {
    auto escape = [](QString s) {
        if (s.contains('"') || s.contains('\t') || s.contains('\n')) {
            s.replace("\"", "\"\"");
            s = '"' + s + '"';
        }
        return s;
    };
    QString text = columns.join('\t') + '\n';
    for (auto row : rows) {
        QStringList c;
        for (auto n : columns)
            c << escape(row[n]);
        text += c.join('\t') + '\n';
    }
    homology::writeText(path, text.toUtf8());
}
static QString section(QString line) {
    if (line.startsWith("Section III"))
        return "genewise";
    if (line.startsWith("Section II"))
        return "spaln";
    if (line.startsWith("Section I."))
        return "raw";
    return {};
}
Catalog catalog(const QString &path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        throw homology::Error("Не удалось прочитать файл последовательностей");
    Catalog c;
    QSet<QString> genes, species, sections;
    QString sec;
    QTextStream in(&f);
    QRegularExpression group("^(?:raw sequence|spaln|genewise) for [^_]+_(\\S+)$");
    while (!in.atEnd()) {
        auto l = in.readLine().trimmed();
        auto s = section(l);
        if (!s.isEmpty()) {
            sec = s;
            sections.insert(s);
        }
        auto m = group.match(l);
        if (m.hasMatch() && sec == "raw" && !genes.contains(m.captured(1))) {
            genes.insert(m.captured(1));
            c.genes << m.captured(1);
        }
        if (l.startsWith('>')) {
            auto sp = l.mid(1).section('|', 0, 0);
            species.insert(sp);
            c.counts[sp]++;
        }
    }
    if (c.genes.isEmpty() || sections.size() != 3)
        throw homology::Error(
            "В TXT должны быть три раздела Section I, II, III и группы генов. Файл "
            "FASTA другого формата не поддерживается.");
    c.species = species.values();
    c.species.sort();
    if (c.species.removeOne("dmel"))
        c.species.prepend("dmel");
    return c;
}
QVector<Record> records(const QString &path, const QSet<QString> &selected) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        throw homology::Error("Не удалось открыть последовательности");
    QTextStream in(&f);
    QString sec, gene;
    Record current;
    QVector<Record> out;
    auto finish = [&] {
        if (!current.header.isEmpty() && (selected.isEmpty() || selected.contains(current.gene)))
            out << current;
        current = {};
    };
    QRegularExpression group("^(?:raw sequence|spaln|genewise) for [^_]+_(\\S+)$");
    while (!in.atEnd()) {
        auto l = in.readLine().trimmed();
        auto s = section(l);
        auto m = group.match(l);
        if (!s.isEmpty()) {
            finish();
            sec = s;
            gene.clear();
        } else if (m.hasMatch()) {
            finish();
            gene = m.captured(1);
        } else if (l.startsWith('>')) {
            finish();
            current = {gene, l.mid(1).section('|', 0, 0), sec, l, l.mid(1).section('|', 1, 1), {}};
        } else if (!current.header.isEmpty() && !l.isEmpty()) {
            l.remove(QRegularExpression("\\s+"));
            current.dna += l.toUpper();
        }
    }
    finish();
    return out;
}
void subset(const QString &source, const QStringList &genes, const QString &output,
            const QSet<QString> &allowedSpecies) {
    QSet<QString> selection(genes.begin(), genes.end());
    auto r = records(source, selection);
    QString text;
    for (auto sec : {QString("raw"), QString("spaln"), QString("genewise")}) {
        text += "Section " +
                (sec == "raw"     ? QString("I. Raw hits")
                 : sec == "spaln" ? QString("II. Spliced align by Spaln")
                                  : QString("III. Spliced align by Genewise")) +
                "\n";
        for (auto gene : genes) {
            text += (sec == "raw" ? QString("raw sequence") : sec) + " for dmel_" + gene + '\n';
            for (auto record : r)
                if (record.gene == gene && record.section == sec &&
                    (allowedSpecies.isEmpty() || allowedSpecies.contains(record.species)))
                    text += record.header + '\n' + record.dna + '\n';
            text += '\n';
        }
    }
    homology::writeText(output, text.toUtf8());
}
QMap<QString, QString> speciesNames() {
    return {{"dmel", "Drosophila melanogaster"},
            {"dsim", "Drosophila simulans"},
            {"dsec", "Drosophila sechellia"},
            {"dyak", "Drosophila yakuba"},
            {"dere", "Drosophila erecta"},
            {"dfic", "Drosophila ficusphila"},
            {"drho", "Drosophila rhopaloa"},
            {"dele", "Drosophila elegans"},
            {"deug", "Drosophila eugracilis"},
            {"dbia", "Drosophila biarmipes"},
            {"dtak", "Drosophila takahashii"},
            {"dkik", "Drosophila kikkawai"},
            {"dana", "Drosophila ananassae"},
            {"dmir", "Drosophila miranda"},
            {"dpse", "Drosophila pseudoobscura"},
            {"dper", "Drosophila persimilis"},
            {"dwil", "Drosophila willistoni"},
            {"sleb", "Scaptodrosophila lebanonensis"},
            {"bdor", "Bactrocera dorsalis"},
            {"blat", "Bactrocera latifrons"},
            {"mdom", "Musca domestica"},
            {"agam", "Anopheles gambiae"},
            {"aaeg", "Aedes aegypti"},
            {"bmor", "Bombyx mori"},
            {"mpha", "Monomorium pharaonis"},
            {"amel", "Apis mellifera"},
            {"aflo", "Apis florea"},
            {"dmag", "Daphnia magna"}};
}
QString displaySpecies(const QString &c) {
    return speciesNames().value(c, c);
}
QString translate(const QString &dna) {
    if (dna.isEmpty() || dna.size() % 3)
        return {};
    QString bases = "TCAG", aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    QMap<QString, QChar> code;
    int i = 0;
    for (auto a : bases)
        for (auto b : bases)
            for (auto c : bases) {
                QString k;
                k += a;
                k += b;
                k += c;
                code[k] = aa[i++];
            }
    QString p;
    for (int j = 0; j < dna.size(); j += 3) {
        auto codon = dna.mid(j, 3);
        if (!code.contains(codon))
            return {};
        p += code[codon];
    }
    if (!p.endsWith('*') || p.left(p.size() - 1).contains('*') || !dna.startsWith("ATG"))
        return {};
    p.chop(1);
    return p;
}
std::shared_ptr<Node> parseTree(const QString &text) {
    int pos = 0;
    auto white = [&] {
        while (pos < text.size() && text[pos].isSpace())
            ++pos;
    };
    std::function<std::shared_ptr<Node>()> parse = [&] {
        white();
        auto n = std::make_shared<Node>();
        if (pos < text.size() && text[pos] == '(') {
            ++pos;
            for (;;) {
                n->children << parse();
                white();
                if (pos >= text.size())
                    throw homology::Error("Незавершённое дерево HAL");
                if (text[pos] == ')') {
                    ++pos;
                    break;
                }
                if (text[pos++] != ',')
                    throw homology::Error("Ошибка формата Newick");
            }
        }
        white();
        while (pos < text.size() && !QString(":,();").contains(text[pos]))
            n->name += text[pos++];
        n->name = n->name.trimmed();
        if (pos < text.size() && text[pos] == ':') {
            ++pos;
            QString length;
            while (pos < text.size() && !QString(",();").contains(text[pos]))
                length += text[pos++];
            bool ok;
            n->length = length.toDouble(&ok);
            if (!ok || !std::isfinite(n->length) || n->length < 0)
                throw homology::Error("Некорректная длина ветви");
        }
        return n;
    };
    auto root = parse();
    white();
    if (pos >= text.size() || text[pos++] != ';')
        throw homology::Error("Дерево должно завершаться точкой с запятой");
    white();
    if (pos != text.size())
        throw homology::Error("Лишние данные после дерева");
    return root;
}
QStringList leafNames(const std::shared_ptr<Node> &n) {
    if (n->children.isEmpty())
        return {n->name};
    QStringList out;
    for (auto c : n->children)
        out << leafNames(c);
    return out;
}
QString className(const QString &c) {
    return QMap<QString, QString>{{"de_novo", "De novo"},
                                  {"duplication", "Duplication"},
                                  {"ancient_diverged", "Ancient / Diverged"},
                                  {"other", "Other"},
                                  {"ambiguous", "Ambiguous"}}
        .value(c, "Ambiguous");
}
QJsonObject classify(const QString &gene, const QJsonObject &e) {
    QJsonArray trace, limitations, alternatives;
    auto rule = [&](QString code, QString text, bool met) {
        trace.append(QJsonObject{{"code", code}, {"text", text}, {"met", met}});
    };
    bool primary = e.value("primary_valid").toBool(),
         complete = e.value("homology_complete").toBool(), orf = e.value("orf_available").toBool(),
         syn = e.value("synteny_available").toBool();
    int paralogs = e.value("strong_paralogs").toInt(),
        distant = e.value("strong_distant_species").toArray().size(),
        homologs = e.value("accepted_cross_species").toArray().size(),
        context = e.value("context_species").toArray().size(),
        broken = e.value("interrupted_species").toArray().size();
    rule("reference", "Подтверждён белковый запрос выбранного референса", primary);
    rule("search", "Все запланированные поиски завершены", complete);
    rule("paralogs", "Сильные совпадения с другими генами референса", paralogs > 0);
    rule("distant", "Сильное сходство не менее чем у двух удалённых видов", distant >= 2);
    rule("synteny", "Есть двусторонняя поддержка геномного окружения", syn && context > 0);
    rule("orf", "Получены сравнительные профили ORF", orf);
    rule("ancestral", "Есть независимая проверка предкового некодирующего состояния",
         e.value("ancestral_noncoding_verified").toBool());
    QString c = "ambiguous",
            reason = "Данных недостаточно для выбора одного механизма происхождения.";
    int score = 15;
    if (!primary) {
        reason = "Для этого кандидата не найден однозначный аннотированный белок "
                 "выбранного основного организма.";
        score = 0;
    } else if (!complete) {
        reason = "Есть незавершённые или ошибочные поиски. Отсутствие результата "
                 "нельзя считать отсутствием гомолога.";
        score = 5;
    } else if (e.value("other_mechanism_verified").toBool()) {
        c = "other";
        reason = "Есть отдельно подтверждённые свидетельства другого механизма "
                 "происхождения.";
        score = 65;
    } else if (!homologs && !paralogs && context >= 2 &&
               e.value("ancestral_noncoding_verified").toBool()) {
        c = "de_novo";
        reason = "Отсутствуют принятые гомологи и паралоги; соответствующий "
                 "предковый некодирующий локус независимо подтверждён.";
        score = 70;
    } else if (paralogs > 0 && distant >= 2) {
        reason = "Одновременно поддерживаются родственные копии и широкое эволюционное "
                 "распространение. Направление и время дупликации не установлены.";
        score = 35;
        alternatives = QJsonArray{"duplication", "ancient_diverged"};
    } else if (paralogs > 0) {
        c = "duplication";
        reason = "Обнаружены сильные кандидаты в паралоги; дупликация "
                 "рассматривается как рабочая гипотеза.";
        score = (syn && context > 0 && orf) ? 60 : 45;
        limitations.append("Общая родственная копия не устанавливает направление и "
                           "время дупликации.");
    } else if (distant >= 2) {
        c = "ancient_diverged";
        reason = "Сильное сходство у нескольких удалённых видов поддерживает более "
                 "широкое происхождение родственной последовательности.";
        score = (syn && context > 0 && orf) ? 60 : 45;
        limitations.append("Ортология и точный возраст гена не установлены.");
    } else if (!homologs && context >= 2 && broken >= 2) {
        reason = "Синтенный контекст и нарушения предсказанных CDS позволяют "
                 "рассматривать De novo, но не исключают потерю кодирования у "
                 "родственников.";
        score = 30;
        alternatives.append("de_novo");
    } else if (homologs > 0) {
        reason = "Есть межвидовое сходство, но оно не определяет механизм "
                 "возникновения гена.";
        score = 30;
    }
    if (e.value("short_query").toBool()) {
        score = qMin(score, 35);
        limitations.append("Короткий запрос: чувствительность и специфичность "
                           "стандартного поиска ограничены.");
    }
    if (e.value("multiple_loci").toBool()) {
        score = qMin(score, 40);
        limitations.append("Есть несколько неоднозначных локусов; один ортолог не выбран.");
    }
    if (!syn)
        limitations.append("Нет применимой оценки синтении.");
    if (!orf)
        limitations.append("Нет пригодного сравнительного профиля ORF.");
    limitations.append("Баллы — эвристическая поддержка рабочей гипотезы, не "
                       "вероятность её истинности.");
    limitations.append("Предсказания Spaln и GeneWise не являются независимыми экспериментами.");
    limitations.append("ORF и синтения не объединяются в один локус без проверки координат.");
    return {{"gene_id", gene},
            {"class", c},
            {"class_label", className(c)},
            {"support_score", score},
            {"confidence_level", score >= 55 ? "Умеренная"
                                 : score > 0 ? "Низкая"
                                             : "Недостаточно данных"},
            {"reason", reason},
            {"rule_trace", trace},
            {"alternatives", alternatives},
            {"limitations", limitations},
            {"evidence", e}};
}
QString findHalBin() {
    const auto bundled = QDir(QCoreApplication::applicationDirPath()).absoluteFilePath("../Tools/hal/bin");
    if (QFileInfo::exists(bundled + "/halStats")) return QDir(bundled).canonicalPath();
    auto path = QStandardPaths::findExecutable("halStats");
    if (!path.isEmpty())
        return QFileInfo(path).absolutePath();
    for (auto p : {QString("/opt/homebrew/Caskroom/miniconda/base/envs/hackathon/bin"),
                   QDir::homePath() + "/miniforge3/envs/hackathon/bin",
                   QDir::homePath() + "/miniconda3/envs/hackathon/bin"})
        if (QFileInfo::exists(p + "/halStats"))
            return p;
    return {};
}
QString findBlastBin(const QString &w) {
    const auto bundled = QDir(QCoreApplication::applicationDirPath()).absoluteFilePath("../Tools/blast/bin");
    if (QFileInfo::exists(bundled + "/blastp")) return QDir(bundled).canonicalPath();
    auto p = w + "/homology/tools/blast/bin";
    if (QFileInfo::exists(p + "/blastp"))
        return p;
    auto t = QStandardPaths::findExecutable("blastp");
    return t.isEmpty() ? QString() : QFileInfo(t).absolutePath();
}
Result loadResult(const QString &dir) {
    auto m = homology::readJson(dir + "/run.json");
    if ((m.value("application") != "LibreGenesis" && m.value("application") != "GeneOrigin") || m.value("status") != "completed")
        throw homology::Error("Каталог не содержит завершённого результата LibreGenesis.");
    for (QString filename : {"summary.json", "summary.tsv", "species_tree.nwk"}) {
        auto hashes = m.value("output_sha256").toObject();
        if (hashes.contains(filename) &&
            homology::sha256(dir + "/" + filename) != hashes.value(filename).toString())
            throw homology::Error("Контрольная сумма файла результата не совпадает: " + filename);
    }
    QFile f(dir + "/summary.json");
    if (!f.open(QIODevice::ReadOnly))
        throw homology::Error("Не найден summary.json");
    QJsonParseError parseError;
    auto data = QJsonDocument::fromJson(f.readAll(), &parseError).array();
    if (parseError.error != QJsonParseError::NoError || data.isEmpty())
        throw homology::Error("Повреждён или пуст summary.json");
    if (!m.value("result_sha256").toString().isEmpty() &&
        homology::sha256(dir + "/summary.json") != m.value("result_sha256").toString())
        throw homology::Error("Контрольная сумма summary.json не совпадает. Результат изменён.");
    QFile tree(dir + "/species_tree.nwk");
    if (!tree.open(QIODevice::ReadOnly))
        throw homology::Error("Не найдено дерево результата");
    return {dir, QString::fromUtf8(tree.readAll()), data, m.value("reference_species").toString(),
            m};
}
} // namespace go
