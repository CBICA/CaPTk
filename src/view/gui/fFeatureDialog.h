#ifndef _fFeatureDialog_h_
#define _fFeatureDialog_h_
//#include "featureParser.h"
#include "FeatureExtraction/src/depends/featureParser.h"

// #include <QtGui>
// NEW CHANGES
#include <QtWidgets>
#include <iostream>
#include "qtreewidget.h"
#include "cbicaUtilities.h"

inline void removeChar(std::string& str, char ch)
{
  str.erase(remove(str.begin(), str.end(), ch), str.end());
}
//TBD add documentation once class is frozen. Also move defs to cpp file
class ParamWidget : public QWidget
{
  Q_OBJECT
public:

  explicit ParamWidget(const std::map< std::string, std::string > field, QWidget *parent = 0) : QWidget(parent)
  {
    QVBoxLayout *mainLayout = new QVBoxLayout;
    QHBoxLayout *hbox = new QHBoxLayout;
    editCombo = NULL;
    editLine = NULL;

    m_widget = getFieldWidget(field);

    m_buttonReset = new QPushButton();
    QCommonStyle* myStyle = new QCommonStyle;
    QIcon openIcon = myStyle->standardIcon(QStyle::SP_BrowserReload);

    m_buttonReset->setIcon(openIcon);
    m_buttonReset->setFixedWidth(30);
    m_buttonReset->setToolTip("Reset");
    m_buttonReset->setEnabled(false);
    QLabel *label = new QLabel(m_widget->accessibleName());
    m_widget->setFixedWidth(100);
    label->setFixedWidth(150);
    hbox->addWidget(label);
    hbox->addWidget(m_widget);
    hbox->addWidget(m_buttonReset);
    label->setWordWrap(true);
    if (editLine)
    {
      connect(editLine, SIGNAL(textChanged(QString)), this, SLOT(enableReset()));
    }
    else if (editCombo)
    {
      connect(editCombo, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(enableReset()));
    }
    connect(m_buttonReset, SIGNAL(clicked()), this, SLOT(resetValue()));
    mainLayout->addLayout(hbox);
    mainLayout->addStretch();
    setLayout(mainLayout);
    enableReset();
  }
  std::string getValue()
  {
    std::string value;
    if (editCombo)
    {
      value = editCombo->currentText().toStdString();
    }
    else if (editLine)
    {
      value = editLine->text().toStdString();
    }
    return value;

  }
  private slots :
  void enableReset()
  {
    m_buttonReset->setEnabled(getValue() != m_defaultValue);
  }
  void resetValue()
  {
    m_buttonReset->setEnabled(false);
    if (editLine)
    {
      editLine->setText(m_defaultValue.c_str());
    }
    else if (editCombo)
    {
      int index = editCombo->findText(m_defaultValue.c_str());
      if (index != -1) { // -1 for not found
        editCombo->setCurrentIndex(index);
      }
    }
  }
private:
  QWidget* getFieldWidget(std::map< std::string, std::string > field)
  {
    std::string name = field["ParamName"];
    std::string value = field["Value"];
    std::string range = field["Range"];
    m_defaultValue = field["Default"];
    std::string type = field["Type"];
    std::transform(type.begin(), type.end(), type.begin(), tolower);
    std::string comment = field["Comments"];
    if (range.find("[") != std::string::npos)
    {
      editCombo = new QComboBox();
    }
    else
    {
      editLine = new QLineEdit();
    }
    removeChar(range, '[');
    removeChar(range, ']');
    removeChar(range, '(');
    removeChar(range, ')');
    auto ranges = cbica::stringSplit(range, ":");
    QWidget* edit = NULL;
    if (editCombo)
    {
      for (auto val : ranges)
      {
        editCombo->addItem(val.c_str());
      }
      int index = editCombo->findText(value.c_str());
      if (index != -1) { // -1 for not found
        editCombo->setCurrentIndex(index);
      }
      edit = editCombo;
    }
    else
    {
      editLine->setText(value.c_str());
      if (type == "int")
      {
        QIntValidator* validater = new QIntValidator(this);
        if (ranges.size() == 2)
        {
          validater->setRange(atoi(ranges[0].c_str()), atoi(ranges[1].c_str()));
        }
        editLine->setValidator(validater);
      }
      else if (type == "real")
      {
        QDoubleValidator* validater = new QDoubleValidator(this);
        if (ranges.size() == 2)
        {
          validater->setRange((double)atof(ranges[0].c_str()), (double)atof(ranges[1].c_str()));
        }
        editLine->setValidator(validater);
      }
      edit = editLine;
    }
    edit->setAccessibleName(name.c_str());
    edit->setToolTip(comment.c_str());
    return edit;
  }

  std::string m_defaultValue;
  QWidget* m_widget;
  QComboBox* editCombo;
  QLineEdit* editLine;
  QPushButton* m_buttonReset;
};
class FeatureDialogTree : public QDialog
{
  Q_OBJECT

public:
  explicit FeatureDialogTree(QWidget *parent = 0) : QDialog(parent)
  {
    treeWidget = new  QTreeWidget();
    treeWidget->setHeaderHidden(true);
    QVBoxLayout *mainLayout = new QVBoxLayout;
    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    m_lblMessage = new QLabel();
    buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    treeWidget->setFixedSize(QSize(400, 400));
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(treeWidget);
    mainLayout->addLayout(buttonsLayout);
    mainLayout->addWidget(m_lblMessage);
    mainLayout->addWidget(buttonBox);
    m_lblMessage->setWordWrap(true);
    m_lblMessage->setVisible(false);
    setLayout(mainLayout);
    setWindowTitle(tr("Feature panel: Advanced parameters"));
    this->setMaximumSize(QSize(400, 400));
  }
  void showMessage(const std::string& message)
  {
    m_lblMessage->setText(message.c_str());
    m_lblMessage->setVisible(true);
  }
  void addFeatureMap(std::map< std::string, std::vector< std::map< std::string, std::string> > > featureMap, std::map< std::string, bool> selectedFeature = std::map< std::string, bool>())
  {
    m_paramWidgets.clear();
    for (auto &features : featureMap)
    {
      if (!selectedFeature.empty())
      {
        if (!selectedFeature[features.first])
        {
          continue;
        }
      }

      if ((features.first == "Intensity") || (features.first == "Laws") || (features.first == "PowerSpectrum"))
      {
        continue;
      }

      QTreeWidgetItem* parentItem = addParent(treeWidget->invisibleRootItem(), features.first);
      auto fields = features.second;

      for (auto field : fields)
      {
        auto tempCheck = field["ParamName"];
        // do not populate for the feature families that do not have an option for parameterization
        if (tempCheck.empty() || (tempCheck == "") || (tempCheck == " "))
        {
          continue;
        }
        QTreeWidgetItem* item = new QTreeWidgetItem(parentItem);
        ParamWidget* wdg = new ParamWidget(field, NULL);
        treeWidget->setItemWidget(item, 0, wdg);
        wdg->setAccessibleName(features.first.c_str());
        m_paramWidgets.push_back(wdg);
      }
    }
    m_featureMap = featureMap;
  }
  
  std::map< std::string, std::vector< std::map< std::string, std::string> > > getFeatureMap()
  {
    return m_featureMap;
  }
  
  public slots:
  void accept()
  {
    updateFeatureMap();
    QDialog::accept();
  }

private:
  void updateFeatureMap()
  {
    std::map< std::string, std::vector< std::string > > updatedVals;
    for (auto p : m_paramWidgets)
    {
      std::string featureName = p->accessibleName().toStdString();
      std::string value = p->getValue();
      updatedVals[featureName].push_back(value);

    }
    for (auto feature : updatedVals)
    {
      for (size_t i = 0; i < feature.second.size(); i++)
      {
        if (i < m_featureMap[feature.first].size())
        {
          m_featureMap[feature.first][i]["Value"] = feature.second[i];
        }
      }
    }
  }
  QTreeWidgetItem* addParent(QTreeWidgetItem *parent, const std::string&  title)
  {
    QTreeWidgetItem*item = new QTreeWidgetItem(parent);
    item->setText(0, tr(title.c_str()));
    item->setExpanded(false);
    return item;
  }
  QTreeWidget *treeWidget;
  QDialogButtonBox *buttonBox;
  QLabel* m_lblMessage;
  std::map< std::string, std::vector< std::map< std::string, std::string> > > m_featureMap;
  std::vector<ParamWidget*> m_paramWidgets;

};
#endif