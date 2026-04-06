#ifndef MONOMIALPARSER_H
#define MONOMIALPARSER_H
#include <QString>
#include <QPair>
class MonomialParser
{
public:
    static QPair<int,int> parse(const QString& mono);
};
#endif
