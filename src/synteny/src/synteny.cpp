#include <synteny/synteny.hpp>

#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSet>
#include <QTemporaryDir>
#include <algorithm>
#include <limits>

namespace synteny {
namespace {

[[noreturn]] void fail(const QString &message) {
    throw homology::Error(message);
}

qint64 number(const QString &text) {
    bool ok = false;
    const auto value = text.toLongLong(&ok);
    if (!ok || value < 0)
        fail("Invalid nonnegative integer: " + text);
    return value;
}

void checkCancel(const homology::Hooks &hooks) {
    if (hooks.cancelled && hooks.cancelled->load())
        fail("cancelled");
}

QString clean(QString text) {
    return text.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
}

void atomicWrite(const QString &path, const QByteArray &bytes) {
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
        fail("Cannot write " + path);
}

using Row = QMap<QString, QString>;

void table(const QString &path, const QStringList &required,
           const std::function<void(const Row &)> &visit) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        fail("Cannot read " + path);
    const auto header = QString::fromUtf8(file.readLine()).trimmed().split('\t');
    for (const auto &column : required) {
        if (!header.contains(column))
            fail(path + ": missing column " + column);
    }
    if (QSet<QString>(header.begin(), header.end()).size() != header.size())
        fail(path + ": duplicate columns");

    qint64 line = 1;
    while (!file.atEnd()) {
        ++line;
        auto text = QString::fromUtf8(file.readLine());
        while (text.endsWith('\n') || text.endsWith('\r'))
            text.chop(1);
        if (text.isEmpty())
            continue;
        const auto fields = text.split('\t');
        if (fields.size() != header.size())
            fail(path + ":" + QString::number(line) + ": column count mismatch");
        Row row;
        for (int i = 0; i < header.size(); ++i)
            row[header[i]] = fields[i];
        visit(row);
    }
}

QVector<qint64> list(QString text, int count) {
    if (text.endsWith(','))
        text.chop(1);
    QVector<qint64> values;
    for (const auto &field : text.split(','))
        values.append(number(field));
    if (values.size() != count)
        fail("PSL block count mismatch");
    return values;
}

QJsonObject flybaseSource(const QJsonObject &manifest, const QString &reference = "dmel") {
    QJsonObject source;
    for (const auto &value : manifest.value("sources").toArray()) {
        const auto entry = value.toObject();
        if (!entry.value("enabled").toBool(true) || entry.value("species") != reference)
            continue;
        if (!source.isEmpty())
            fail("Ambiguous dmel FlyBase sources in sources_manifest");
        source = entry;
    }
    if (source.isEmpty() || (reference == "dmel" && !source.value("assembly").toString().startsWith("dmel_r6.36_")))
        fail("Unverified dmel r6.36 source assembly");
    return source;
}

QJsonObject selectionProvenance(const QJsonObject &config, const QDir &base,
                                const QString &annotation, const QString &selectionPath,
                                const QJsonObject &selection, const QJsonObject &source,
                                QJsonObject &inputHashes) {
    const auto selectedHash = inputHashes.value("annotation_tsv").toString();
    if (!QFile::exists(selectionPath)) {
        return {{"status", "selection_not_provided"},
                {"selected_annotation_sha256", selectedHash}};
    }

    QString metadataPath;
    if (config.contains("source_metadata")) {
        const auto configured = config.value("source_metadata").toString();
        if (configured.isEmpty())
            fail("Empty synteny.source_metadata");
        metadataPath = base.absoluteFilePath(configured);
    } else {
        const auto id = source.value("id").toString();
        if (!QRegularExpression("^[A-Za-z0-9_][A-Za-z0-9_.-]*$").match(id).hasMatch())
            fail("Missing/invalid FlyBase source id; specify synteny.source_metadata");
        metadataPath = QFileInfo(annotation).dir().absoluteFilePath("../" + id + "/source.json");
    }
    metadataPath = QDir::cleanPath(metadataPath);
    const auto metadata = homology::readJson(metadataPath);
    const auto specification = metadata.value("specification").toObject();
    for (const auto &key : {"species", "provider", "assembly", "id"}) {
        const auto expected = source.value(key).toString();
        if (expected.isEmpty() || specification.value(key).toString() != expected)
            fail("Source metadata specification mismatch: " + QString(key));
    }

    QMap<QString, QString> records;
    const QRegularExpression sha256Pattern("^[0-9a-fA-F]{64}$");
    for (const auto &value : metadata.value("files").toArray()) {
        const auto record = value.toObject();
        const auto name = record.value("name").toString();
        if (name != "annotation.tsv" && name != "proteins.faa")
            continue;
        const auto digest = record.value("sha256").toString();
        if (records.contains(name) || !sha256Pattern.match(digest).hasMatch())
            fail("Invalid/duplicate source metadata hash: " + name);
        records[name] = digest.toLower();
    }
    for (const auto &pair : {qMakePair(QString("annotation.tsv"), QString("annotation_sha256")),
                             qMakePair(QString("proteins.faa"), QString("proteome_sha256"))}) {
        const auto expected = selection.value(pair.second).toString();
        if (!sha256Pattern.match(expected).hasMatch() || !records.contains(pair.first) ||
            records.value(pair.first) != expected.toLower())
            fail("Selection source metadata hash mismatch: " + pair.first);
    }

    const auto metadataHash = homology::sha256(metadataPath);
    inputHashes["source_metadata"] = metadataHash;
    // extractQueries hashes the full inputs, not the selected output files. Metadata
    // can attest their recorded hashes, but cannot verify the absent full inputs.
    return {{"status", "provenance_metadata_only"},
            {"source_metadata_path", metadataPath},
            {"source_metadata_sha256", metadataHash},
            {"source_specification", specification},
            {"source_annotation_sha256", records.value("annotation.tsv")},
            {"source_proteome_sha256", records.value("proteins.faa")},
            {"selected_annotation_sha256", selectedHash},
            {"full_source_files_verified", false},
            {"limitations", "Selection source hashes matched recorded metadata only; "
                            "full source annotation/proteome contents were not verified"}};
}

struct Gene {
    QString id, protein;
    Location loc;
    qint64 length = 0;
};

struct Feature {
    QString gene, kind, sequence;
    Interval interval;
};

struct Hsp {
    QString gene, species, assembly, sequence, quality;
    QChar strand;
    Interval interval;
};

struct Group {
    QString sequence;
    QChar orientation;
    QVector<Block> cds, left, right;
};

QVector<Interval> sources(const QVector<Block> &blocks) {
    QVector<Interval> intervals;
    for (const auto &block : blocks)
        intervals.append(block.source);
    return intervals;
}

Interval bounds(const QVector<Block> &blocks) {
    Interval interval{std::numeric_limits<qint64>::max(), 0};
    for (const auto &block : blocks) {
        interval.start = std::min(interval.start, block.target.start);
        interval.end = std::max(interval.end, block.target.end);
    }
    return interval;
}

qint64 cov(const QVector<Block> &blocks) {
    return coveredLength(sources(blocks));
}

bool discordant(QVector<Block> blocks) {
    std::sort(blocks.begin(), blocks.end(), [](const auto &a, const auto &b) {
        return a.source.start < b.source.start;
    });
    for (int i = 1; i < blocks.size(); ++i) {
        const auto &a = blocks[i - 1];
        const auto &b = blocks[i];
        if (a.source.end > b.source.start)
            return true;
        if (a.orientation == '+' ? a.target.end > b.target.start : b.target.end > a.target.start)
            return true;
    }
    return false;
}

QMap<QString, Group> groups(const QVector<Block> &blocks,
                           const QMap<QString, Feature> &features) {
    QMap<QString, Group> result;
    for (const auto &block : blocks) {
        if (!features.contains(block.name))
            fail("Unknown PSL BED name: " + block.name);
        const auto feature = features[block.name];
        if (block.sourceSequence != feature.sequence ||
            block.source.start < feature.interval.start || block.source.end > feature.interval.end)
            fail("PSL source outside input feature");
        const auto key = feature.gene + '\t' + block.sequence + '\t' + block.orientation;
        auto &group = result[key];
        group.sequence = block.sequence;
        group.orientation = block.orientation;
        if (feature.kind == "CDS")
            group.cds.append(block);
        else if (feature.kind == "left")
            group.left.append(block);
        else
            group.right.append(block);
    }
    return result;
}

bool supported(const Group &group, qint64 minimum) {
    return cov(group.left) >= minimum && cov(group.right) >= minimum &&
           flankOrder(bounds(group.left), bounds(group.right), group.orientation);
}

QStringList csv(const QString &text) {
    QStringList fields;
    QString field;
    bool quoted = false;
    for (int i = 0; i < text.size(); ++i) {
        const auto ch = text[i];
        if (ch == '"') {
            if (quoted && i + 1 < text.size() && text[i + 1] == '"') {
                field += '"';
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (ch == ',' && !quoted) {
            fields << field;
            field.clear();
        } else {
            field += ch;
        }
    }
    if (quoted)
        fail("Unterminated CSV quote");
    fields << field;
    return fields;
}

QString joinRow(const QStringList &fields) {
    QStringList result;
    for (const auto &field : fields)
        result << clean(field);
    return result.join('\t') + '\n';
}

} // namespace

qint64 coveredLength(QVector<Interval> intervals) {
    for (const auto interval : intervals) {
        if (interval.start < 0 || interval.end < interval.start)
            fail("Invalid interval");
    }
    std::sort(intervals.begin(), intervals.end(), [](auto a, auto b) {
        return a.start < b.start;
    });
    qint64 total = 0;
    qint64 end = 0;
    for (const auto interval : intervals) {
        total += std::max<qint64>(0, interval.end - std::max(end, interval.start));
        end = std::max(end, interval.end);
    }
    return total;
}

Location parseLocation(const QString &text) {
    Location location;
    const auto colon = text.indexOf(':');
    if (colon <= 0)
        fail("Invalid location: " + text);
    location.sequence = text.left(colon);
    QString expression = text.mid(colon + 1);
    expression.remove(QRegularExpression("\\s+"));
    if (expression.startsWith("complement(") && expression.endsWith(')')) {
        location.strand = '-';
        expression = expression.mid(11, expression.size() - 12);
    }
    if (expression.startsWith("join(") && expression.endsWith(')'))
        expression = expression.mid(5, expression.size() - 6);

    const QRegularExpression range("^([0-9]+)\\.\\.([0-9]+)$");
    for (const auto &part : expression.split(',')) {
        const auto match = range.match(part);
        if (!match.hasMatch())
            fail("Invalid location: " + text);
        const auto start = number(match.captured(1));
        const auto end = number(match.captured(2));
        if (start < 1 || end < start)
            fail("Invalid location: " + text);
        location.blocks.append({start - 1, end});
    }
    // FlyBase complement(join(...)) can list exons in descending genomic order.
    std::sort(location.blocks.begin(), location.blocks.end(), [](auto a, auto b) {
        return a.start < b.start;
    });
    for (int i = 1; i < location.blocks.size(); ++i) {
        if (location.blocks[i].start < location.blocks[i - 1].end)
            fail("Overlapping CDS exons: " + text);
    }
    return location;
}

QVector<Block> parsePsl(const QString &line) {
    auto fields = line.trimmed().split(QRegularExpression("\\s+"));
    QString name;
    if (fields.size() == 22)
        name = fields.takeFirst();
    else if (fields.size() != 21)
        fail("PSL requires 21 fields plus optional leading BED name");
    for (int i = 0; i < 8; ++i)
        number(fields[i]);

    const auto strand = fields[8];
    if (strand != "+" && strand != "++" && strand != "+-")
        fail("PSL expected source +, target +/-");
    const auto querySize = number(fields[10]);
    const auto targetSize = number(fields[14]);
    const auto queryStart = number(fields[11]);
    const auto queryEnd = number(fields[12]);
    const auto targetStart = number(fields[15]);
    const auto targetEnd = number(fields[16]);
    if (queryStart > queryEnd || queryEnd > querySize || targetStart > targetEnd ||
        targetEnd > targetSize)
        fail("Invalid PSL bounds");
    const auto count = number(fields[17]);
    if (count < 1 || count > 10000000)
        fail("Invalid PSL block count");
    const auto sizes = list(fields[18], count);
    const auto queryStarts = list(fields[19], count);
    const auto targetStarts = list(fields[20], count);
    const QChar orientation = strand.endsWith('-') ? '-' : '+';
    QVector<Block> blocks;
    for (int i = 0; i < count; ++i) {
        const auto size = sizes[i];
        if (size <= 0 || size > querySize || size > targetSize ||
            queryStarts[i] > querySize - size || targetStarts[i] > targetSize - size)
            fail("Invalid PSL block coordinates");
        const auto start = orientation == '-' ? targetSize - targetStarts[i] - size : targetStarts[i];
        if (queryStarts[i] < queryStart || queryStarts[i] + size > queryEnd ||
            start < targetStart || start + size > targetEnd)
            fail("PSL block outside alignment bounds");
        blocks.append({name, fields[13], orientation,
                       {queryStarts[i], queryStarts[i] + size}, {start, start + size}, fields[9]});
    }
    return blocks;
}

bool flankOrder(Interval left, Interval right, QChar orientation) {
    return orientation == '+' ? left.end <= right.start : orientation == '-' && right.end <= left.start;
}

bool compatible(const QString &assembly, const QString &otherAssembly,
                const QString &sequence, const QString &otherSequence,
                QChar strand, QChar otherStrand) {
    return !assembly.isEmpty() && assembly == otherAssembly && sequence == otherSequence &&
           strand == otherStrand;
}

QJsonObject run(const QString &configPath, const homology::Hooks &hooks) {
    checkCancel(hooks);
    const auto root = homology::readJson(configPath);
    const auto config = root.value("synteny").toObject();
    if (config.isEmpty())
        fail("Missing synteny config");
    const QDir base(QFileInfo(configPath).absolutePath());
    const auto path = [&](const QString &key) {
        const auto value = config.value(key).toString();
        if (value.isEmpty())
            fail("Missing synteny." + key);
        return base.absoluteFilePath(value);
    };
    const auto outputDir = path("output_dir");
    if (QFile::exists(outputDir + "/run.json"))
        fail("Existing synteny run.json: " + outputDir);
    if (QDir(outputDir).exists() &&
        !QDir(outputDir).entryList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot).isEmpty())
        fail("Nonempty output directory: " + outputDir);

    const auto homologyDir = path("homology_dir");
    const auto annotation = path("annotation_tsv");
    const auto query = path("query_fasta");
    const auto queryMap = path("query_map");
    const auto hal = path("hal_file");
    const auto homologyRun = homology::readJson(homologyDir + "/run.json");
    if (homologyRun.value("status") != "completed")
        fail("Homology run not completed");

    QJsonObject inputHashes;
    const auto hash = [&](const QString &key, const QString &file) {
        const auto digest = homology::sha256(file);
        inputHashes[key] = digest;
        return digest;
    };
    if (hash("query_fasta", query) != homologyRun.value("query_sha256").toString() ||
        hash("query_map", queryMap) != homologyRun.value("query_map_sha256").toString())
        fail("Homology query hashes mismatch");
    const auto outputHashes = homologyRun.value("output_sha256").toObject();
    for (const auto &file : {"homology_hits.tsv", "homology_summary.tsv", "candidate_loci.tsv",
                             "search_status.tsv"}) {
        if (hash(file, homologyDir + '/' + file) != outputHashes.value(file).toString())
            fail(QString("Homology hash mismatch: ") + file);
    }
    hash("annotation_tsv", annotation);
    hash("homology_run", homologyDir + "/run.json");
    hash("sources_manifest", path("sources_manifest"));
    hash("assembly_csv", path("assembly_csv"));
    hash("pipeline_config", configPath);
    const auto sourceManifest = homology::readJson(path("sources_manifest"));
    const auto referenceSpecies = config.value("reference_species").toString("dmel");
    const auto sourceSpecification = flybaseSource(sourceManifest, referenceSpecies);
    const auto selectionPath = QFileInfo(annotation).dir().filePath("selection.json");
    QJsonObject selection;
    if (QFile::exists(selectionPath)) {
        selection = homology::readJson(selectionPath);
        hash("selection", selectionPath);
    }
    const auto provenance = selectionProvenance(config, base, annotation, selectionPath,
                                                selection, sourceSpecification, inputHashes);

    QMap<QString, QString> proteinGene;
    table(queryMap, {"protein_id", "gene_id"}, [&](const Row &row) {
        if (proteinGene.contains(row["protein_id"]))
            fail("Duplicate query map protein");
        proteinGene[row["protein_id"]] = row["gene_id"];
    });
    QMap<QString, QString> proteins;
    for (const auto &sequence : homology::readFasta(query, true)) {
        if (proteins.contains(sequence.id))
            fail("Duplicate query FASTA ID");
        proteins[sequence.id] = sequence.sequence;
    }
    QMap<QString, Gene> all;
    table(annotation, {"protein_id", "gene_id", "transcript_id", "species", "location", "length"},
          [&](const Row &row) {
        const auto id = row["gene_id"];
        const auto protein = row["protein_id"];
        if (id.isEmpty() || protein.isEmpty() || row["transcript_id"].isEmpty() ||
            row["species"] != referenceSpecies || all.contains(id))
            fail("Invalid/duplicate annotation gene " + id);
        if (proteinGene.value(protein) != id || !proteins.contains(protein) ||
            number(row["length"]) != proteins[protein].size())
            fail("Selection mismatch for " + id);
        Gene gene;
        gene.id = id;
        gene.protein = protein;
        try {
            gene.loc = parseLocation(row["location"]);
        } catch (const std::exception &error) {
            fail(id + ": " + error.what());
        }
        gene.length = coveredLength(gene.loc.blocks);
        all[id] = gene;
    });
    if (all.size() != proteinGene.size() || proteins.size() != all.size() ||
        all.size() != homologyRun.value("query_count").toInt())
        fail("Annotation/query count mismatch");
    if (QFile::exists(selectionPath)) {
        QSet<QString> selected;
        for (const auto &value : selection.value("genes").toArray()) {
            const auto id = value.toString();
            if (!all.contains(id) || selected.contains(id))
                fail("Invalid selection gene " + id);
            selected.insert(id);
        }
        if (selected.size() != all.size())
            fail("Selection gene set mismatch");
    }
    QMap<QString, Gene> genes = all;
    if (config.contains("gene_ids")) {
        genes.clear();
        for (const auto &value : config.value("gene_ids").toArray()) {
            const auto id = value.toString();
            if (!all.contains(id) || genes.contains(id))
                fail("Invalid requested gene " + id);
            genes[id] = all[id];
        }
        if (genes.isEmpty())
            fail("Empty gene_ids");
    }

    QMap<QString, QString> assemblies;
    for (const auto &value : homologyRun.value("config").toObject().value("targets").toArray()) {
        const auto target = value.toObject();
        const auto species = target.value("species").toString();
        if (species == referenceSpecies || species == "aaeg" || !target.value("enabled").toBool(true))
            continue;
        const auto assembly = target.value("assembly").toString();
        if (assemblies.contains(species) && assemblies[species] != assembly)
            fail("Multiple assemblies for " + species);
        assemblies[species] = assembly;
    }
    if (assemblies.isEmpty())
        fail("No comparison species configured");

    QMap<QString, QString> assemblyFull;
    QFile assemblyCsv(path("assembly_csv"));
    if (!assemblyCsv.open(QIODevice::ReadOnly | QIODevice::Text))
        fail("Cannot read assembly CSV");
    if (csv(QString::fromUtf8(assemblyCsv.readLine()).trimmed()) != QStringList{"Species", "Assembly"})
        fail("Invalid assembly CSV header");
    while (!assemblyCsv.atEnd()) {
        const auto fields = csv(QString::fromUtf8(assemblyCsv.readLine()).trimmed());
        if (fields.size() != 2)
            fail("Invalid assembly CSV row");
        if (assemblyFull.contains(fields[1]))
            fail("Duplicate assembly CSV mapping");
        assemblyFull[fields[1]] = QString(fields[0]).replace(' ', '_');
        if (fields[0] == "Drosophila melanogaster" && fields[1].contains("r6.36"))
            assemblyFull["dmel_r6.36_FB2020_05"] = "Drosophila_melanogaster";
    }
    QMap<QString, QString> full;
    QSet<QString> verified;
    for (const auto &value : sourceManifest.value("sources").toArray()) {
        const auto source = value.toObject();
        if (!source.value("enabled").toBool(true))
            continue;
        const auto species = source.value("species").toString();
        const auto assembly = source.value("assembly").toString();
        if (assemblies.value(species) == assembly && assemblyFull.contains(assembly)) {
            full[species] = assemblyFull[assembly];
            verified.insert(species);
        }
    }
    QString source = assemblyFull.value(sourceSpecification.value("assembly").toString());
    if (source.isEmpty()) fail("Reference assembly not present in assembly CSV");

    QMap<QString, QVector<Hsp>> hsps;
    table(homologyDir + "/homology_hits.tsv",
          {"gene_id", "target_species", "assembly", "search_type", "target_sequence_id",
           "target_start", "target_end", "strand", "passes_filters", "quality"},
          [&](const Row &row) {
        checkCancel(hooks);
        if (!genes.contains(row["gene_id"]) || !assemblies.contains(row["target_species"]) ||
            !row["search_type"].contains("tblastn"))
            return;
        if (row["passes_filters"] != "true" && row["quality"] != "significant_partial_similarity")
            return;
        const auto start = number(row["target_start"]);
        const auto end = number(row["target_end"]);
        if (start < 1 || end < 1 || (row["strand"] != "+" && row["strand"] != "-"))
            fail("Invalid HSP coordinates");
        hsps[row["gene_id"] + '\t' + row["target_species"]].append(
            {row["gene_id"], row["target_species"], row["assembly"], row["target_sequence_id"],
             row["quality"], row["strand"][0], {std::min(start, end) - 1, std::max(start, end)}});
    });
    QMap<QString, QString> geneGroups;
    table(homologyDir + "/homology_summary.tsv",
          {"gene_id", "protein_hit_species", "genome_hit_species"}, [&](const Row &row) {
        const bool protein = !row["protein_hit_species"].isEmpty();
        const bool genome = !row["genome_hit_species"].isEmpty();
        geneGroups[row["gene_id"]] = !protein && !genome ? "no_accepted"
                                  : !protein && genome ? "tblastn_only" : "others";
    });
    table(homologyDir + "/search_status.tsv", {"gene_id", "status"}, [&](const Row &row) {
        if (!row["status"].startsWith("completed"))
            fail("Incomplete homology search for " + row["gene_id"]);
    });

    const int primary = config.value("primary_flank").toInt(2000);
    const int fallback = config.value("fallback_flank").toInt(10000);
    const int minimum = config.value("min_flank_bases").toInt(100);
    const int timeout = config.value("timeout_seconds").toInt(600);
    if (primary < 1 || fallback < primary || minimum < 1 || timeout < 1)
        fail("Invalid flank/timeout configuration");
    if (!QDir().mkpath(outputDir))
        fail("Cannot create output");
    QTemporaryDir temporary;
    if (!temporary.isValid())
        fail("Cannot create temporary directory");

    QJsonArray commands, errors, scope, unsupported;
    for (auto it = assemblies.begin(); it != assemblies.end(); ++it)
        scope.append(QJsonObject{{"species", it.key()}, {"assembly", it.value()}});
    QJsonObject manifest{
        {"resolved_parameters", QJsonObject{{"primary_flank", primary}, {"fallback_flank", fallback},
                                            {"min_flank_bases", minimum}, {"timeout_seconds", timeout}}},
        {"schema_version", 1},
        {"module_version", "0.1.0"},
        {"status", "running"},
        {"config", config},
        {"input_hashes", inputHashes},
        {"selection_provenance", provenance},
        {"coordinate_system", "BED0 half-open; source location is CDS, not full gene"},
        {"flank_definition", "genomic left/right independent of CDS strand; input BED strand always +"},
        {"scope", scope},
        {"excluded_aaeg", true},
        {"gene_count", genes.size()},
        {"started_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
        {"limitations", QJsonArray{
            "No origin classification or automatic orthology/duplication calls",
            "Context intervals are not ancestral CDS",
            "HSPs are not gene boundaries",
            "Same-contig dispersed mappings retained conservatively; no unique ortholog assignment",
            "HAL assembly verification relies on supplied assembly CSV version"}}
    };
    const auto save = [&]() {
        manifest["commands"] = commands;
        manifest["errors"] = errors;
        manifest["unsupported"] = unsupported;
        atomicWrite(outputDir + "/run.json", QJsonDocument(manifest).toJson());
    };
    save();
    const auto exec = [&](const QString &tool, const QStringList &arguments) {
        checkCancel(hooks);
        const auto executable = QDir(path("hal_bin")).filePath(tool);
        QJsonArray args;
        for (const auto &argument : arguments)
            args.append(argument);
        commands.append(QJsonObject{{"executable", executable}, {"arguments", args}});
        return homology::process(executable, arguments, timeout, hooks);
    };

    QByteArray results =
        "gene_id\tspecies\tassembly\thal_genome\tstatus\tmessage\tsource_strand\tcds_bases\tcds_covered\t"
        "primary_flank\tused_flank\tprimary_left_covered\tprimary_right_covered\tleft_covered\t"
        "right_covered\tleft_available\tright_available\tcontext_groups\tmapping_groups\thomology_hsps\t"
        "target_orientations\tmapping_complexity\n";
    QByteArray loci =
        "gene_id\tspecies\tassembly\tcontig\tstart_0based\tend_0based_exclusive\tstrand\tprovenance\t"
        "selection_justification\tquality_limitations\n";
    QByteArray evidence =
        "gene_id\tspecies\thsp_assembly\tcontig\tstart_0based\tend_0based_exclusive\tstrand\tquality\t"
        "comparison_status\tdirect_cds_overlap_bases\tcontext_overlap_bases\n";
    QMap<QString, int> counts, mappedSpecies;
    try {
        QJsonObject versions;
        for (const auto &tool : {"halStats", "halLiftover"}) {
            QJsonObject version{
                {"executable_sha256", homology::sha256(QDir(path("hal_bin")).filePath(tool))}};
            QString banner;
            bool successfulProbe = false;
            try {
                banner = QString::fromUtf8(exec(tool, {"--version"})).trimmed();
                successfulProbe = true;
            } catch (const std::exception &error) {
                checkCancel(hooks);
                banner = QString::fromUtf8(error.what());
            }
            const QRegularExpression pattern(
                "\\b" + QRegularExpression::escape(tool) +
                    "\\s*(?:version\\s*[:=]?\\s*|v\\s*)?"
                    "([0-9]+(?:\\.[0-9]+)+(?:[-+][A-Za-z0-9_.-]+)?)\\b",
                QRegularExpression::CaseInsensitiveOption);
            const auto match = pattern.match(banner);
            if (match.hasMatch()) {
                version["version"] = match.captured(1);
                version["version_detection"] = successfulProbe
                    ? "successful_probe_banner" : "nonzero_probe_banner";
            } else {
                version["version_detection"] = "unavailable";
                version["version_probe_error"] = clean(
                    (successfulProbe ? "Version not found in --version output: " : "") + banner).left(512);
            }
            versions[tool] = version;
        }
        manifest["versions"] = versions;
        const auto genomes = QString::fromUtf8(exec("halStats", {hal, "--genomes"}))
                                 .split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        if (!genomes.contains(source))
            fail("Source genome unavailable in HAL");
        QMap<QString, qint64> chromSizes;
        const auto sizes = QString::fromUtf8(exec("halStats", {hal, "--chromSizes", source}));
        for (const auto &line : sizes.split('\n', Qt::SkipEmptyParts)) {
            const auto fields = line.simplified().split(' ');
            if (fields.size() != 2)
                fail("Invalid halStats chromSizes");
            chromSizes[fields[0]] = number(fields[1]);
        }
        for (const auto &gene : genes) {
            if (!chromSizes.contains(gene.loc.sequence) ||
                gene.loc.blocks.last().end > chromSizes[gene.loc.sequence])
                fail("Source CDS outside HAL chromosome: " + gene.id);
        }

        for (auto species = assemblies.begin(); species != assemblies.end(); ++species) {
            checkCancel(hooks);
            if (hooks.progress)
                hooks.progress("SYNTENY " + species.key());
            QString status = "completed";
            QString message;
            QString fallbackMessage;
            bool fallbackSucceeded = false;
            QMap<QString, Feature> features;
            QVector<Block> primaryBlocks, fallbackBlocks;
            QSet<QString> needsFallback;
            QMap<QString, Group> primaryGroups, fallbackGroups;
            const auto batch = [&](int flank, const QSet<QString> &subset, bool cds) {
                const auto bed = temporary.filePath("input.bed");
                const auto psl = temporary.filePath("output.psl");
                QFile::remove(psl);
                QByteArray bytes;
                int index = features.size();
                for (const auto &gene : genes) {
                    if (!subset.contains(gene.id))
                        continue;
                    const auto add = [&](const QString &kind, Interval interval) {
                        if (interval.end <= interval.start)
                            return;
                        const auto name = "f" + QString::number(index++);
                        features[name] = {gene.id, kind, gene.loc.sequence, interval};
                        bytes += joinRow({gene.loc.sequence, QString::number(interval.start),
                                          QString::number(interval.end), name, "0", "+"}).toUtf8();
                    };
                    if (cds) {
                        for (const auto block : gene.loc.blocks)
                            add("CDS", block);
                    }
                    add("left", {std::max<qint64>(0, gene.loc.blocks.first().start - flank),
                                 gene.loc.blocks.first().start});
                    add("right", {gene.loc.blocks.last().end,
                                  std::min<qint64>(chromSizes[gene.loc.sequence],
                                                   gene.loc.blocks.last().end + flank)});
                }
                atomicWrite(bed, bytes);
                exec("halLiftover", {hal, source, bed, full[species.key()], psl, "--outPSLWithName"});
                QFile file(psl);
                if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
                    fail("Missing PSL output");
                QVector<Block> blocks;
                while (!file.atEnd()) {
                    checkCancel(hooks);
                    const auto line = QString::fromUtf8(file.readLine()).trimmed();
                    if (!line.isEmpty())
                        blocks += parsePsl(line);
                }
                return blocks;
            };

            if (!verified.contains(species.key())) {
                status = "assembly_unverified";
            } else if (!genomes.contains(full[species.key()])) {
                status = "unavailable_in_hal";
            } else {
                try {
                    QSet<QString> subset;
                    for (const auto &gene : genes)
                        subset.insert(gene.id);
                    primaryBlocks = batch(primary, subset, true);
                    primaryGroups = groups(primaryBlocks, features);
                    for (const auto &gene : genes) {
                        bool both = false;
                        for (auto it = primaryGroups.begin(); it != primaryGroups.end(); ++it) {
                            if (it.key().section('\t', 0, 0) == gene.id && supported(it.value(), minimum))
                                both = true;
                        }
                        if (!both)
                            needsFallback.insert(gene.id);
                    }
                } catch (const std::exception &error) {
                    checkCancel(hooks);
                    status = "liftover_failed";
                    message = QString::fromUtf8(error.what());
                    errors.append(QJsonObject{{"species", species.key()}, {"stage", "primary"},
                                              {"error", message}});
                    primaryBlocks.clear();
                    primaryGroups.clear();
                }
                if (status == "completed" && !needsFallback.isEmpty() && fallback > primary) {
                    try {
                        fallbackBlocks = batch(fallback, needsFallback, false);
                        fallbackGroups = groups(fallbackBlocks, features);
                        fallbackSucceeded = true;
                    } catch (const std::exception &error) {
                        checkCancel(hooks);
                        const auto detail = QString::fromUtf8(error.what());
                        fallbackMessage = "fallback_liftover_failed;using_primary_window: " + detail;
                        QJsonArray affectedGenes;
                        for (const auto &gene : genes) {
                            if (needsFallback.contains(gene.id))
                                affectedGenes.append(gene.id);
                        }
                        errors.append(QJsonObject{{"species", species.key()}, {"stage", "fallback"},
                                                  {"gene_ids", affectedGenes}, {"error", detail}});
                        fallbackBlocks.clear();
                        fallbackGroups.clear();
                    }
                }
            }
            if (status != "completed")
                unsupported.append(QJsonObject{{"species", species.key()}, {"status", status}});

            for (const auto &gene : genes) {
                const bool used = needsFallback.contains(gene.id) && fallbackSucceeded;
                const auto rowMessage = needsFallback.contains(gene.id) && !fallbackMessage.isEmpty()
                    ? fallbackMessage : message;
                QMap<QString, Group> current;
                for (auto it = primaryGroups.begin(); it != primaryGroups.end(); ++it) {
                    if (it.key().section('\t', 0, 0) != gene.id)
                        continue;
                    auto group = it.value();
                    if (used) {
                        group.left.clear();
                        group.right.clear();
                    }
                    current[it.key()] = group;
                }
                if (used) {
                    for (auto it = fallbackGroups.begin(); it != fallbackGroups.end(); ++it) {
                        if (it.key().section('\t', 0, 0) != gene.id)
                            continue;
                        auto &group = current[it.key()];
                        group.sequence = it->sequence;
                        group.orientation = it->orientation;
                        group.left = it->left;
                        group.right = it->right;
                    }
                }
                QVector<Interval> cdsSources, primaryLeft, primaryRight, usedLeft, usedRight;
                for (auto it = primaryGroups.begin(); it != primaryGroups.end(); ++it) {
                    if (it.key().section('\t', 0, 0) == gene.id) {
                        primaryLeft += sources(it->left);
                        primaryRight += sources(it->right);
                    }
                }
                struct Region {
                    QString sequence;
                    QChar strand;
                    Interval interval;
                    QString provenance;
                };
                QVector<Region> regions;
                int contexts = 0;
                QStringList orientations;
                bool complex = current.size() > 1;
                for (const auto &group : current) {
                    if (!orientations.contains(QString(group.orientation)))
                        orientations << QString(group.orientation);
                    QVector<Block> ordered = group.left;
                    ordered += group.cds;
                    ordered += group.right;
                    const auto span = bounds(ordered);
                    const auto sourceSpan = gene.loc.blocks.last().end - gene.loc.blocks.first().start;
                    if (discordant(ordered) ||
                        span.end - span.start > std::max<qint64>(
                            100000, 10 * (sourceSpan + 2 * (used ? fallback : primary))))
                        complex = true;
                    cdsSources += sources(group.cds);
                    usedLeft += sources(group.left);
                    usedRight += sources(group.right);
                    const QChar strand = gene.loc.strand == group.orientation ? '+' : '-';
                    for (const auto &block : group.cds)
                        regions.append({group.sequence, strand, block.target, "direct_CDS"});
                    if (supported(group, minimum)) {
                        ++contexts;
                        const auto left = bounds(group.left);
                        const auto right = bounds(group.right);
                        const Interval interval = group.orientation == '+'
                            ? Interval{left.end, right.start} : Interval{right.end, left.start};
                        if (interval.end > interval.start)
                            regions.append({group.sequence, strand, interval, "context_only"});
                    }
                }

                const bool hasContext = std::any_of(regions.begin(), regions.end(), [](const auto &region) {
                    return region.provenance == "context_only";
                });
                const auto hits = hsps.value(gene.id + '\t' + species.key());
                for (const auto &hit : hits) {
                    QVector<Interval> direct, context;
                    const bool assemblyOK = verified.contains(species.key()) &&
                                            species.value() == hit.assembly;
                    for (const auto &region : regions) {
                        if (!assemblyOK || !compatible(species.value(), hit.assembly, region.sequence,
                                                       hit.sequence, region.strand, hit.strand))
                            continue;
                        const Interval overlap{std::max(region.interval.start, hit.interval.start),
                                               std::min(region.interval.end, hit.interval.end)};
                        if (overlap.end > overlap.start)
                            (region.provenance == "direct_CDS" ? direct : context).append(overlap);
                    }
                    const auto directCoverage = coveredLength(direct);
                    const auto contextCoverage = coveredLength(context);
                    QString comparison;
                    QString comparisonLimitation;
                    if (!assemblyOK) {
                        comparison = "assembly_unverified_or_mismatch";
                        comparisonLimitation = "assembly_coordinates_not_comparable";
                    } else if (status != "completed") {
                        comparison = status == "unavailable_in_hal" ? "context_unavailable_in_hal"
                            : "mapping_unavailable_" + status;
                        comparisonLimitation = "mapping_status=" + status + ";context_not_assessed";
                    } else if (regions.isEmpty()) {
                        comparison = "context_not_established";
                        comparisonLimitation = "context_not_established;location_relative_to_context_unknown";
                    } else if (directCoverage) {
                        comparison = "direct_CDS_overlap";
                    } else if (contextCoverage) {
                        comparison = "context_overlap";
                    } else if (!hasContext) {
                        comparison = "outside_mapped_CDS_without_context";
                        comparisonLimitation = "outside_mapped_CDS;context_not_established;"
                                               "location_relative_to_context_unknown";
                    } else {
                        comparison = "outside_selected_regions";
                        comparisonLimitation = "outside_selected_regions;not_proof_of_nonorthology";
                    }
                    evidence += joinRow({gene.id, species.key(), hit.assembly, hit.sequence,
                                         QString::number(hit.interval.start), QString::number(hit.interval.end),
                                         QString(hit.strand), hit.quality, comparison,
                                         QString::number(directCoverage), QString::number(contextCoverage)}).toUtf8();
                    if (directCoverage == 0 && contextCoverage == 0) {
                        loci += joinRow({gene.id, species.key(), hit.assembly, hit.sequence,
                                         QString::number(hit.interval.start), QString::number(hit.interval.end),
                                         QString(hit.strand), "homology_only", comparison,
                                         "HSP_not_gene_boundary;alternative_evidence_not_duplication;" +
                                             comparisonLimitation}).toUtf8();
                    }
                }
                for (const auto &region : regions) {
                    loci += joinRow({gene.id, species.key(), species.value(), region.sequence,
                                     QString::number(region.interval.start), QString::number(region.interval.end),
                                     QString(region.strand), region.provenance,
                                     region.provenance == "direct_CDS" ? "mapped_CDS_block"
                                         : "inferred_context_between_collinear_flanks",
                                     complex ? "candidate_only;orthology_unresolved;complex_discordant_or_dispersed_mappings"
                                         : "candidate_only;orthology_unresolved;blocks_may_be_disconnected"}).toUtf8();
                }
                const auto cdsCoverage = coveredLength(cdsSources);
                QString rowStatus = status;
                if (status == "completed") {
                    if (cdsCoverage > 0)
                        rowStatus = complex ? "direct_CDS_complex" : "direct_CDS_candidate";
                    else if (contexts)
                        rowStatus = "inferred_context";
                    else if (!current.isEmpty())
                        rowStatus = "one_sided_or_discordant_context";
                    else
                        rowStatus = "no_alignment";
                }
                ++counts[rowStatus];
                if (cdsCoverage || contexts)
                    ++mappedSpecies[gene.id];
                const auto flank = used ? fallback : primary;
                results += joinRow({
                    gene.id, species.key(), species.value(), full.value(species.key()), rowStatus, rowMessage,
                    QString(gene.loc.strand), QString::number(gene.length), QString::number(cdsCoverage),
                    QString::number(primary), QString::number(flank),
                    QString::number(coveredLength(primaryLeft)), QString::number(coveredLength(primaryRight)),
                    QString::number(coveredLength(usedLeft)), QString::number(coveredLength(usedRight)),
                    QString::number(std::min<qint64>(flank, gene.loc.blocks.first().start)),
                    QString::number(std::min<qint64>(flank, chromSizes[gene.loc.sequence] - gene.loc.blocks.last().end)),
                    QString::number(contexts), QString::number(current.size()), QString::number(hits.size()),
                    orientations.join(','), complex ? "multiple_discordant_or_dispersed_mappings"
                        : "not_detected_not_proof_of_uniqueness"
                }).toUtf8();
            }
            save();
        }
        manifest["status"] = errors.isEmpty() ? "completed" : "partial";
    } catch (const std::exception &error) {
        manifest["status"] = hooks.cancelled && hooks.cancelled->load() ? "cancelled" : "failed";
        manifest["fatal_error"] = QString::fromUtf8(error.what());
    }

    QByteArray summary = "gene_id\tprotein_id\thomology_group\tspecies_with_CDS_or_context\t"
                         "scope_species\tlimitations\n";
    for (const auto &gene : genes) {
        summary += joinRow({gene.id, gene.protein, geneGroups.value(gene.id, "unavailable"),
                            QString::number(mappedSpecies[gene.id]), QString::number(assemblies.size()),
                            "no_origin_classification;orthology_unresolved"}).toUtf8();
    }
    QJsonObject hashes;
    const QVector<QPair<QString, QByteArray>> outputs{
        {"synteny_results.tsv", results}, {"synteny_loci.tsv", loci},
        {"gene_summary.tsv", summary}, {"synteny_evidence.tsv", evidence}
    };
    for (const auto &output : outputs) {
        atomicWrite(outputDir + '/' + output.first, output.second);
        hashes[output.first] = homology::sha256(outputDir + '/' + output.first);
    }
    manifest["output_hashes"] = hashes;
    QJsonObject statusCounts;
    for (auto it = counts.begin(); it != counts.end(); ++it)
        statusCounts[it.key()] = it.value();
    manifest["counts"] = statusCounts;
    manifest["finished_utc"] = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
    manifest["hal_file"] = QJsonObject{
        {"path", hal}, {"size", double(QFileInfo(hal).size())},
        {"mtime_utc", QFileInfo(hal).lastModified().toUTC().toString(Qt::ISODate)},
        {"sha256_status", "not_hashed_large_input"}
    };
    save();
    return manifest;
}

} // namespace synteny
