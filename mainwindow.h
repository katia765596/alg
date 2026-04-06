#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QGraphicsScene>
#include <QVector>
#include <QPair>
class QLineEdit;
class QPushButton;
class QGraphicsView;
class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
private slots:
    void onUpdateClicked();
private:
    void setupUI();
    void drawPoints(const QVector<QPair<int,int>>& goodPoints,
                    const QVector<QPair<int,int>>& badPoints);
    QLineEdit *lineEdit;
    QPushButton *btn;
    QGraphicsView *view;
    QGraphicsScene *scene;
    static constexpr int gridSize = 20;
    static constexpr int cellSize = 30;
};
#endif
