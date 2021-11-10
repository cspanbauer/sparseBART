/*
 *  sparseBART: sparse Bayesian Additive Regression Trees
 *  Copyright (C) 2021 Charles Spanbauer
 *
 *  This file is part of sparseBART.
 *
 *  sparseBART is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; version 3 of the License, or
 *  (at your option) any later version.
 *
 *  sparseBART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with sparseBART; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-3
 */

#include <BART3.h>

#ifdef build_with_DPM
#include "DPM.h"
#include "DPMneal7.h"
#include "DPMneal8.h"
#include "DPMNoGa.h"
#endif

#include "cEXPVALUE.h"
#include "cgbart.h"
#include "cpwbart.h"
#include "chotdeck.h"
#include "mc_cores_openmp.h"
#include "RcppEigen.h"
