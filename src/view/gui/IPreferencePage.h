#ifndef IPREFERENCEPAGE_H
#define IPREFERENCEPAGE_H

#include <QObject>
class IPreferencePage
{
protected:
    IPreferencePage();

public slots:
    virtual void OnOkay() = 0;
    virtual void OnCancel() = 0;
	virtual void Restore() = 0;

};

#endif // IPREFERENCEPAGE_H

Q_DECLARE_INTERFACE(IPreferencePage, "IPreferencePage");
