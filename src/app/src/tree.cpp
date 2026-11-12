#include "tree.hpp"
#include <QGraphicsScene>
#include <QGraphicsTextItem>
#include <QMouseEvent>
#include <QPainter>
#include <QSet>
#include <QResizeEvent>
#include <QWheelEvent>
#include <functional>
EvidenceTree::EvidenceTree(QWidget *p) : QGraphicsView(p), scene_(new QGraphicsScene(this)) {
    setScene(scene_);
    setRenderHint(QPainter::Antialiasing);
    setDragMode(ScrollHandDrag);
    setTransformationAnchor(NoAnchor);
    setBackgroundBrush(QColor("#FAFCFB"));
    setFrameShape(QFrame::NoFrame);
    setMinimumHeight(260);
}
double EvidenceTree::minimumScale() const {
    const auto rect = scene_->sceneRect();
    if (rect.isEmpty()) return 1;
    return qMin(viewport()->width() / (rect.width() + 4),
                viewport()->height() / (rect.height() + 4));
}
void EvidenceTree::fitTree() {
    if (scene_->items().isEmpty()) return;
    autoFit_ = true;
    resetTransform();
    const double fitted = minimumScale();
    scale(fitted, fitted);
    centerOn(scene_->sceneRect().center());
}
void EvidenceTree::changeZoom(double factor) {
    if (scene_->items().isEmpty()) return;
    const auto rect = scene_->sceneRect();
    const auto focus = mapToScene(viewport()->rect().center());
    const double current = transform().m11();
    const double next = qBound(minimumScale(), current * factor, qMax(4.0, minimumScale()));
    autoFit_ = qAbs(next - minimumScale()) < .001;
    resetTransform();
    scale(next, next);
    // Keep the visible area inside the tree, centering axes that fit entirely.
    const double halfWidth = viewport()->width() / (2.0 * next);
    const double halfHeight = viewport()->height() / (2.0 * next);
    const double x = 2 * halfWidth >= rect.width() ? rect.center().x()
        : qBound(rect.left() + halfWidth, focus.x(), rect.right() - halfWidth);
    const double y = 2 * halfHeight >= rect.height() ? rect.center().y()
        : qBound(rect.top() + halfHeight, focus.y(), rect.bottom() - halfHeight);
    centerOn(x, y);
}
void EvidenceTree::zoomIn() { changeZoom(1.2); }
void EvidenceTree::zoomOut() { changeZoom(1 / 1.2); }
void EvidenceTree::resizeEvent(QResizeEvent *event) {
    QGraphicsView::resizeEvent(event);
    if (autoFit_) fitTree();
    else changeZoom(1);
}
void EvidenceTree::wheelEvent(QWheelEvent *e) {
    if (e->angleDelta().y() > 0)
        zoomIn();
    else
        zoomOut();
    e->accept();
}
void EvidenceTree::mousePressEvent(QMouseEvent *e) {
    auto item = itemAt(e->pos());
    while (item) {
        auto sp = item->data(0).toString();
        if (!sp.isEmpty()) {
            if (highlight_) {
                scene_->removeItem(highlight_);
                delete highlight_;
                highlight_ = nullptr;
            }
            double y = item->sceneBoundingRect().center().y();
            for (auto *node : scene_->items())
                if (node->data(0).toString() == sp && dynamic_cast<QGraphicsEllipseItem *>(node)) {
                    y = node->sceneBoundingRect().center().y();
                    break;
                }
            highlight_ = scene_->addRect(QRectF(480, y - 23, scene_->sceneRect().right() - 480, 45),
                                         Qt::NoPen, QColor("#E2EFE7"));
            highlight_->setZValue(-1);
            emit speciesSelected(sp);
            break;
        }
        item = item->parentItem();
    }
    QGraphicsView::mousePressEvent(e);
}
void EvidenceTree::setResult(const QString &newick, const QJsonObject &ev, const QString &gene,
                             const QString &reference) {
    scene_->clear();
    highlight_ = nullptr;
    primary_ = reference;
    if (newick.isEmpty())
        return;
    auto root = go::parseTree(newick);
    auto names = go::speciesNames();
    QMap<QString, QString> codes;
    for (auto it = names.begin(); it != names.end(); ++it)
        codes[QString(it.value()).replace(' ', '_')] = it.key();
    auto species = ev.value("species").toObject();
    QSet<QString> present;
    std::function<std::shared_ptr<go::Node>(std::shared_ptr<go::Node>)> prune = [&](auto n) {
        auto x = std::make_shared<go::Node>(*n);
        x->children.clear();
        if (n->children.isEmpty()) {
            if (species.contains(codes.value(n->name))) {
                present.insert(codes.value(n->name));
                return x;
            }
            return std::shared_ptr<go::Node>();
        }
        for (auto c : n->children) {
            auto q = prune(c);
            if (q)
                x->children << q;
        }
        if (x->children.size() == 1)
            return x->children.first();
        return x->children.isEmpty() ? std::shared_ptr<go::Node>() : x;
    };
    root = prune(root);
    if (!root)
        return;
    int maxDepth = 0;
    std::function<void(std::shared_ptr<go::Node>, int)> depth = [&](auto n, int d) {
        maxDepth = qMax(maxDepth, d);
        for (auto c : n->children)
            depth(c, d + 1);
    };
    depth(root, 0);
    double row = 45;
    int leaf = 0;
    double step = 430. / qMax(1, maxDepth);
    QPen branch(QColor("#C6D6D1"), 1.8);
    std::function<QPointF(std::shared_ptr<go::Node>, int)> draw = [&](auto n, int d) {
        double x = 32 + d * step;
        if (n->children.isEmpty()) {
            double y = 55 + leaf++ * row;
            auto code = codes.value(n->name);
            auto data = species.value(code).toObject();
            QString orf = data.value("orf").toString(), hom = data.value("homology").toString();
            QColor color("#A9B9B3");
            if (hom == "принятое совпадение")
                color = QColor("#4382B4");
            if (hom == "частичное сходство")
                color = QColor("#8E79AC");
            if (orf.startsWith("preserved_both"))
                color = QColor("#168D77");
            if (orf == "predicted_interruption_both" || orf == "discordant_methods" ||
                orf == "conflicting_alignment")
                color = QColor("#D8963D");
            if (code == reference)
                color = QColor("#163F3B");
            double end = 470;
            scene_->addLine(x, y, end, y, branch);
            auto dot = scene_->addEllipse(end - 5, y - 5, 10, 10, QPen(color, 2), color);
            dot->setData(0, code);
            if (code == reference) {
                auto *ring = scene_->addEllipse(end - 9, y - 9, 18, 18, QPen(QColor("#D8AE61"), 2));
                ring->setData(0, code);
                ring->setZValue(-0.1);
            }
            auto text = scene_->addText(
                go::displaySpecies(code) + (code == reference ? "  • основной" : ""),
                QFont("Arial", 10, code == reference ? QFont::DemiBold : QFont::Normal));
            text->setPos(488, y - 20);
            text->setDefaultTextColor(QColor("#24443D"));
            text->setData(0, code);
            QString target = code == reference ? ev.value("reference_gene_id").toString()
                                               : data.value("target_label").toString();
            if (target.isEmpty())
                target = code == reference ? ev.value("reference_gene_id").toString()
                                           : "локус не установлен";
            auto small = scene_->addText(target + "  ·  " + gene, QFont("Arial", 8));
            small->setPos(488, y + 1);
            small->setDefaultTextColor(QColor("#80938C"));
            small->setData(0, code);
            QString tip = go::displaySpecies(code) + "\nГомология: " + hom +
                          "\nСинтения: " + data.value("synteny").toString() + "\nORF: " + orf;
            dot->setToolTip(tip);
            text->setToolTip(tip);
            small->setToolTip(tip);
            return QPointF(x, y);
        }
        QVector<QPointF> points;
        for (auto c : n->children)
            points << draw(c, d + 1);
        double y = (points.first().y() + points.last().y()) / 2;
        scene_->addLine(x, points.first().y(), x, points.last().y(), branch);
        for (auto p : points)
            scene_->addLine(x, p.y(), p.x(), p.y(), branch);
        auto dot =
            scene_->addEllipse(x - 2.5, y - 2.5, 5, 5, QPen(QColor("#B2C6BD")), QColor("#B2C6BD"));
        dot->setToolTip("Узел дерева HAL. Предковый ген и кодирующее состояние "
                        "здесь не реконструированы.");
        return QPointF(x, y);
    };
    draw(root, 0);
    auto title = scene_->addText("ДЕРЕВО ВИДОВ · HAL", QFont("Arial", 9, QFont::DemiBold));
    title->setPos(30, 7);
    title->setDefaultTextColor(QColor("#66827A"));
    QStringList missing;
    for (auto key : species.keys())
        if (!present.contains(key))
            missing << go::displaySpecies(key);
    double bottom = 55 + leaf * row;
    if (!missing.isEmpty()) {
        auto m = scene_->addText("Вне дерева HAL: " + missing.join(", ") +
                                     "\nТопологическое положение не дорисовывается. "
                                     "Данные этих видов доступны в таблице.",
                                 QFont("Arial", 9));
        m->setTextWidth(760);
        m->setPos(30, bottom);
        m->setDefaultTextColor(QColor("#8C8068"));
    }
    scene_->setSceneRect(scene_->itemsBoundingRect().adjusted(-20, -12, 25, 20));
    fitTree();
}
bool EvidenceTree::savePng(const QString &path) {
    auto rect = scene_->sceneRect();
    double factor = qMin(2.0, 6000.0 / qMax(rect.width(), rect.height()));
    QImage image((rect.size() * factor).toSize(), QImage::Format_ARGB32_Premultiplied);
    image.fill(Qt::white);
    QPainter p(&image);
    p.setRenderHint(QPainter::Antialiasing);
    scene_->render(&p, QRectF(QPointF(0, 0), image.size()), rect);
    p.end();
    return image.save(path);
}
