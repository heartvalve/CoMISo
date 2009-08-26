//=============================================================================
//
//  CLASS MiSolverDialog
//
//=============================================================================


#ifndef ACG_MISOLVERDIALOG_HH
#define ACG_MISOLVERDIALOG_HH


//== INCLUDES =================================================================

// qmake users have to include
#ifdef QT4_FOUND 
#include "ui_QtMISolverDialogBase.hh"

// ACGMake users have to include
// #include "QtMISolverDialogBase.hh"


#include <Solver/MISolver.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG
{

//== CLASS DEFINITION =========================================================



/** \class MISolverDialog MISolverDialog.hh <ACG/.../MISolverDialog.hh>

    Brief Description.

    A more elaborate description follows.
*/
class MISolverDialog
         : public QDialog, public Ui::QtMISolverDialogBase
{
   Q_OBJECT
public:

   /// Default constructor
   MISolverDialog( MISolver&   _misolver,
                   QWidget*    _parent = 0 ):
         QDialog( _parent ),
         Ui::QtMISolverDialogBase(),
         misolver_( _misolver )
   {
      setupUi( this );
      get_parameters();

      connect( okPB, SIGNAL( clicked() ),     this, SLOT( slotOkButton() ) );
      connect( cancelPB, SIGNAL( clicked() ), this, SLOT( slotCancelButton() ) );
   }

   /// Destructor
   ~MISolverDialog() {}

   void get_parameters();
   void set_parameters();


public slots:
   virtual void slotOkButton();
   virtual void slotCancelButton();

private:

   MISolver& misolver_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_MISOLVERDIALOG_HH defined
//=============================================================================

#endif
