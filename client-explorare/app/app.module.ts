import { NgModule }       from '@angular/core';
import { BrowserModule }  from '@angular/platform-browser';
import { FormsModule }    from '@angular/forms';
import { HttpModule }     from '@angular/http';

/*** config ***/
import { APP_CONFIG, AppConfig } from './app.config';

/*** components ***/
import { AppComponent }             from './app.component';
import { PhenoMarksComponent }      from './pheno-marks.component';
import { SimilarDiseasesComponent } from './similar-diseases.component';
import { PhenoPanelComponent }      from './pheno-panel.component';

/*** services ***/
import { PhenotypeService } from './phenotype.service';

/*** routing ***/
import { AppRoutingModule }     from './app-routing.module';

/*** utils ***/
import { Ng2AutoCompleteModule } from 'ng2-auto-complete';

/*** clipboard ***/
import { ClipboardModule } from 'ngx-clipboard';

@NgModule({
  imports: [
    BrowserModule,
    FormsModule,
    HttpModule,
    AppRoutingModule,
    Ng2AutoCompleteModule,
    ClipboardModule
  ],
  declarations: [
    AppComponent,
    PhenoMarksComponent,
    SimilarDiseasesComponent,
    PhenoPanelComponent
  ],
  providers: [ PhenotypeService,
              { provide: APP_CONFIG, useValue: AppConfig } ],
  bootstrap: [ AppComponent ]
})
export class AppModule { }


/*
Copyright 2016 Google Inc. All Rights Reserved.
Use of this source code is governed by an MIT-style license that
can be found in the LICENSE file at http://angular.io/license
*/
