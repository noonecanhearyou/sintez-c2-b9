#pragma once
#include "core.hpp"
#include <QGraphicsView>
class EvidenceTree final : public QGraphicsView {
    Q_OBJECT
  public:
    explicit EvidenceTree(QWidget *parent = nullptr);
    void setResult(const QString &newick, const QJsonObject &evidence, const QString &gene,
                   const QString &reference);
    void fitTree();
    void zoomIn();
    void zoomOut();
    bool savePng(const QString &path);
  signals:
    void speciesSelected(QString code);

  protected:
    void wheelEvent(QWheelEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void resizeEvent(QResizeEvent *) override;

  private:
    QGraphicsScene *scene_;
    double minimumScale() const;
    void changeZoom(double factor);
    bool autoFit_ = true;
    QString primary_;
    QGraphicsRectItem *highlight_ = nullptr;
};
