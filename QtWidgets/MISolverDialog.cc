/*===========================================================================*\
 *                                                                           *
 *                              CoMISo                                       *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                      www.graphics.rwth-aachen.de                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is a part of CoMISo.                                          *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 1    $                                                       *
 *   $Author: zimmer $                                                       *
 *   $Date: 2009-08-05 16:35:37 +0200 (Wed, 05 Aug 2009) $                   *
 *                                                                           *
\*===========================================================================*/


//=============================================================================
//
//  CLASS MISolverDialog - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

#include "MISolverDialog.hh"

#ifdef QT4_FOUND
#include <QtGui>

//== NAMESPACES ===============================================================

namespace ACG {

//== IMPLEMENTATION ========================================================== 


//-----------------------------------------------------------------------------


void 
MISolverDialog::
get_parameters()
{
//     QDoubleSpinBox *localErrorDSB;
//     QSpinBox *localItersSB;
//     QDoubleSpinBox *cgErrorDSB;
//     QSpinBox *cgItersSB;
//     QDoubleSpinBox *taucsErrorDSB;
//     QCheckBox *finalTaucsCB;
//     QCheckBox *initialTaucsCB;
//     QSpinBox *infoSB;
//     QCheckBox *directRoundingCB;

  initialFullCB  ->setChecked( misolver_.get_inital_full());
  finalFullCB    ->setChecked( misolver_.get_final_full());
  directRoundingCB->setChecked( misolver_.get_direct_rounding());
  noRoundingCB    ->setChecked( misolver_.get_no_rounding());

  localItersSB ->setValue( misolver_.get_local_iters());
  localErrorDSB->setValue( log(misolver_.get_local_error())/log(10.0f));

  cgItersSB ->setValue( misolver_.get_cg_iters());
  cgErrorDSB->setValue( log(misolver_.get_cg_error())/log(10.0f));
  
  fullErrorDSB->setValue( log(misolver_.get_full_error())/log(10.0f));

  infoSB->setValue( misolver_.get_noise());
  solverStatsCheckBox->setChecked( misolver_.get_stats( ));
}


//-----------------------------------------------------------------------------


void 
MISolverDialog::
set_parameters()
{
  misolver_.set_inital_full   ( initialFullCB  ->isChecked() );
  misolver_.set_final_full    ( finalFullCB    ->isChecked() );
  misolver_.set_direct_rounding( directRoundingCB->isChecked());
  misolver_.set_no_rounding( noRoundingCB->isChecked());

  misolver_.set_local_iters( localItersSB ->value());
  misolver_.set_local_error( pow(10, localErrorDSB->value()));

  misolver_.set_cg_iters( cgItersSB ->value());
  misolver_.set_cg_error( pow(10, cgErrorDSB->value()));

  misolver_.set_full_error( pow(10, fullErrorDSB->value()));

  misolver_.set_noise( infoSB->value());
  misolver_.set_stats( solverStatsCheckBox->isChecked());
}


//-----------------------------------------------------------------------------


void 
MISolverDialog::
slotOkButton()
{
  set_parameters();
  close();
}


//-----------------------------------------------------------------------------


void
MISolverDialog::
slotCancelButton()
{
  close();
}



//=============================================================================
} // namespace ACG
//=============================================================================

#endif
