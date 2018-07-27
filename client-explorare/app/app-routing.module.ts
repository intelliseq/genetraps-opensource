import { NgModule }             from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

import { PhenoMarksComponent }   from './pheno-marks.component';
import { SimilarDiseasesComponent }   from './similar-diseases.component';
import { PhenoPanelComponent }   from './pheno-panel.component';

const routes: Routes = [
  { path: '', redirectTo: '/pheno-marks', pathMatch: 'full' },
  { path: 'pheno-marks',  component: PhenoMarksComponent },
  { path: 'similar-diseases',  component: SimilarDiseasesComponent },
  { path: 'pheno-panel',     component: PhenoPanelComponent }
];

@NgModule({
  imports: [ RouterModule.forRoot(routes) ],
  exports: [ RouterModule ]
})
export class AppRoutingModule {}


/*
Copyright 2016 Google Inc. All Rights Reserved.
Use of this source code is governed by an MIT-style license that
can be found in the LICENSE file at http://angular.io/license
*/
