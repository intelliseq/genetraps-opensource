"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
var __metadata = (this && this.__metadata) || function (k, v) {
    if (typeof Reflect === "object" && typeof Reflect.metadata === "function") return Reflect.metadata(k, v);
};
var core_1 = require('@angular/core');
var platform_browser_1 = require('@angular/platform-browser');
var forms_1 = require('@angular/forms');
var http_1 = require('@angular/http');
/*** config ***/
var app_config_1 = require('./app.config');
/*** components ***/
var app_component_1 = require('./app.component');
var pheno_marks_component_1 = require('./pheno-marks.component');
var similar_diseases_component_1 = require('./similar-diseases.component');
var pheno_panel_component_1 = require('./pheno-panel.component');
/*** services ***/
var phenotype_service_1 = require('./phenotype.service');
/*** routing ***/
var app_routing_module_1 = require('./app-routing.module');
/*** utils ***/
var ng2_auto_complete_1 = require('ng2-auto-complete');
/*** clipboard ***/
var ngx_clipboard_1 = require('ngx-clipboard');
var AppModule = (function () {
    function AppModule() {
    }
    AppModule = __decorate([
        core_1.NgModule({
            imports: [
                platform_browser_1.BrowserModule,
                forms_1.FormsModule,
                http_1.HttpModule,
                app_routing_module_1.AppRoutingModule,
                ng2_auto_complete_1.Ng2AutoCompleteModule,
                ngx_clipboard_1.ClipboardModule
            ],
            declarations: [
                app_component_1.AppComponent,
                pheno_marks_component_1.PhenoMarksComponent,
                similar_diseases_component_1.SimilarDiseasesComponent,
                pheno_panel_component_1.PhenoPanelComponent
            ],
            providers: [phenotype_service_1.PhenotypeService,
                { provide: app_config_1.APP_CONFIG, useValue: app_config_1.AppConfig }],
            bootstrap: [app_component_1.AppComponent]
        }), 
        __metadata('design:paramtypes', [])
    ], AppModule);
    return AppModule;
}());
exports.AppModule = AppModule;
/*
Copyright 2016 Google Inc. All Rights Reserved.
Use of this source code is governed by an MIT-style license that
can be found in the LICENSE file at http://angular.io/license
*/
//# sourceMappingURL=app.module.js.map