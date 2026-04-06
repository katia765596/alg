#include "mainwindow.h"
#include "monomialparser.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QGraphicsTextItem>
#include <QPen>
#include <QBrush>
#include <QMessageBox>
#include <QDebug>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    setupUI();
}
MainWindow::~MainWindow() {}
void MainWindow::setupUI()
{
    QWidget *central = new QWidget(this);
    setCentralWidget(central);
    QVBoxLayout *mainLayout = new QVBoxLayout(central);
    QHBoxLayout *inputLayout = new QHBoxLayout();
    QLabel *label = new QLabel("Образующие (через запятую):");
    lineEdit = new QLineEdit();
    lineEdit->setText("x^6, x^2y^3, xy^7");
    btn = new QPushButton("Обновить");
    inputLayout->addWidget(label);
    inputLayout->addWidget(lineEdit);
    inputLayout->addWidget(btn);
    mainLayout->addLayout(inputLayout);
    scene = new QGraphicsScene(this);
    view = new QGraphicsView(scene);
    view->setRenderHint(QPainter::Antialiasing);
    view->setFixedSize(cellSize * gridSize + 260, cellSize * gridSize + 70);
    view->setStyleSheet("background-color: #fff0f5;");
    mainLayout->addWidget(view);
    connect(btn, &QPushButton::clicked, this, &MainWindow::onUpdateClicked);
    onUpdateClicked();
}
void MainWindow::onUpdateClicked()
{
    QString text = lineEdit->text();
    QStringList parts = text.split(',', Qt::SkipEmptyParts);
    QVector<QPair<int,int>> gens;
    bool ok = true;
    for (const QString& part : parts) {
        QString trimmed = part.trimmed();
        if (trimmed.isEmpty()) continue;
        QPair<int,int> p = MonomialParser::parse(trimmed);
        qDebug() << "Parsed:" << trimmed << "->" << p.first << p.second;
        if (p.first < 0 || p.second < 0) {
            QMessageBox::warning(this, "Ошибка ввода",
                                 QString("Не удалось разобрать моном: %1").arg(trimmed));
            ok = false;
            break;
        }
        gens.append(p);
    }
    if (!ok) return;
    if (gens.isEmpty()) {
        QMessageBox::information(this, "Нет образующих",
                                 "Задайте хотя бы один моном (например, x^2y)");
        scene->clear();
        return;
    }
    qDebug() << "Generators count:" << gens.size();
    QVector<QPair<int,int>> goodPoints, badPoints;
    for (int m = 0; m < gridSize; ++m) {
        for (int n = 0; n < gridSize; ++n) {
            bool inIdeal = false;
            for (const auto& g : gens) {
                if (m >= g.first && n >= g.second) {
                    inIdeal = true;
                    break;
                }
            }
            if (inIdeal)
                goodPoints.append(qMakePair(m,n));
            else
                badPoints.append(qMakePair(m,n));
        }
    }
    qDebug() << "Ideal points:" << goodPoints.size();
    qDebug() << "Remainder points:" << badPoints.size();
    drawPoints(goodPoints, badPoints);
}
void MainWindow::drawPoints(const QVector<QPair<int,int>>& goodPoints,
                            const QVector<QPair<int,int>>& badPoints)
{
    scene->clear();
    QBrush brushGood(QColor(255, 200, 210));
    QBrush brushBad(QColor(255, 140, 170));
    QPen penGrid(QColor(255, 220, 230));
    QPen penAxis(QColor(150, 0, 80));
    QFont fontNum("Arial", 9);
    fontNum.setBold(true);
    for (int x = 0; x <= gridSize; ++x) {
        int px = x * cellSize;
        scene->addLine(px, 0, px, gridSize * cellSize, penGrid);
        if (x < gridSize) {
            QGraphicsTextItem *t = scene->addText(QString::number(x), fontNum);
            t->setDefaultTextColor(penAxis.color());
            t->setPos(px + 4, -22);
        }
    }
    for (int y = 0; y <= gridSize; ++y) {
        int py = y * cellSize;
        scene->addLine(0, py, gridSize * cellSize, py, penGrid);
        if (y < gridSize) {
            QGraphicsTextItem *t = scene->addText(QString::number(y), fontNum);
            t->setDefaultTextColor(penAxis.color());
            t->setPos(-30, py - 6);
        }
    }
    QFont fontLabel("Arial", 11);
    fontLabel.setBold(true);
    QGraphicsTextItem *xLabel = scene->addText("m (показатель x)", fontLabel);
    xLabel->setDefaultTextColor(penAxis.color());
    xLabel->setPos(gridSize * cellSize / 2 - 50, gridSize * cellSize + 12);
    QGraphicsTextItem *yLabel = scene->addText("n (показатель y)", fontLabel);
    yLabel->setDefaultTextColor(penAxis.color());
    yLabel->setPos(-60, gridSize * cellSize / 2);
    yLabel->setRotation(-90);
    for (const auto& p : goodPoints) {
        int x = p.first * cellSize;
        int y = p.second * cellSize;
        scene->addRect(x, y, cellSize-1, cellSize-1, QPen(Qt::NoPen), brushGood);
    }
    for (const auto& p : badPoints) {
        int x = p.first * cellSize;
        int y = p.second * cellSize;
        scene->addRect(x, y, cellSize-1, cellSize-1, QPen(Qt::NoPen), brushBad);
    }
    QFont fontLeg("Arial", 10);
    fontLeg.setBold(true);
    int legX = gridSize * cellSize + 20;
    int legY = 20;
    scene->addRect(legX, legY, 20, 20, QPen(Qt::NoPen), brushGood);
    QGraphicsTextItem *legGood = scene->addText("Мономы в I", fontLeg);
    legGood->setDefaultTextColor(penAxis.color());
    legGood->setPos(legX + 25, legY + 2);
    scene->addRect(legX, legY + 35, 20, 20, QPen(Qt::NoPen), brushBad);
    QGraphicsTextItem *legBad = scene->addText("Мономы остатка", fontLeg);
    legBad->setDefaultTextColor(penAxis.color());
    legBad->setPos(legX + 25, legY + 37);
}
