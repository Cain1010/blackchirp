#ifndef CHIRPCONFIGWIDGET_H
#define CHIRPCONFIGWIDGET_H

#include <QWidget>

#include "chirpconfig.h"
#include "chirptablemodel.h"

class QSpinBox;

namespace Ui {
class ChirpConfigWidget;
}

class ChirpConfigWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ChirpConfigWidget(QWidget *parent = 0);
    ~ChirpConfigWidget();

    ChirpConfig getChirpConfig();
    QSpinBox *numChirpsBox() const;

public slots:
    void initializeFromSettings();
    void enableEditing(bool enabled);
    void setButtonStates();

    void addSegment();
    void addEmptySegment();
    void insertSegment();
    void insertEmptySegment();
    void moveSegments(int direction);
    void removeSegments();
    void clear();
    void load();

    void updateChirpPlot();

signals:
    void chirpConfigChanged();


private:
    Ui::ChirpConfigWidget *ui;
    ChirpTableModel *p_ctm;
    ChirpConfig d_currentChirpConfig;

    bool isSelectionContiguous(QModelIndexList l);
    void clearList(bool replot=true);
    void updateChirpConfig();


};

#endif // CHIRPCONFIGWIDGET_H
