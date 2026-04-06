#include "monomialparser.h"
#include <QRegularExpression>
#include <QDebug>
QPair<int,int> MonomialParser::parse(const QString& mono)
{
    QString s = mono.simplified().remove(' ');
    if (s.isEmpty()) return qMakePair(-1, -1);
    QRegularExpression re(R"(^(?:x(?:\^(\d+))?)?(?:y(?:\^(\d+))?)?$)",
                          QRegularExpression::CaseInsensitiveOption);
    QRegularExpressionMatch match = re.match(s);
    if (!match.hasMatch()) {
        qDebug() << "No match for:" << s;
        return qMakePair(-1, -1);
    }
    int a = 0, b = 0;
    QString xExp = match.captured(1);
    QString yExp = match.captured(2);
    if (s.contains('x', Qt::CaseInsensitive)) {
        a = xExp.isEmpty() ? 1 : xExp.toInt();
    }
    if (s.contains('y', Qt::CaseInsensitive)) {
        b = yExp.isEmpty() ? 1 : yExp.toInt();
    }
    if (!s.contains('x', Qt::CaseInsensitive) && !s.contains('y', Qt::CaseInsensitive))
        return qMakePair(-1, -1);
    return qMakePair(a, b);
}
