/*===========================================================================*\
 *                                                                           *
 *                               CoMISo                                      *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of CoMISo.                                             *
 *                                                                           *
 *  CoMISo is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  CoMISo is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with CoMISo.  If not, see <http://www.gnu.org/licenses/>.          *
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
