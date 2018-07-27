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
var phenotype_service_1 = require('./phenotype.service');
var AppComponent = (function () {
    function AppComponent(phenotypeService) {
        var _this = this;
        phenotypeService.phenotypes.subscribe(function (value) { _this.phenotypeCount = value.length; });
        phenotypeService.diseases.subscribe(function (value) { _this.diseaseCount = value.length; });
        phenotypeService.genes.subscribe(function (value) { _this.geneCount = value.length; });
    }
    AppComponent = __decorate([
        core_1.Component({
            selector: 'my-app',
            template: "\n<div class=\"container\">\n  <div class=\"row py-1\" data-spy=\"affix\">\n    <div class=\"col-sm-4\">\n      <h1>explorare</h1>\n    </div>\n    <div class=\"col-sm-2\">\n      Phenotypes <span class=\"badge badge-pill badge-success\">{{phenotypeCount}}</span>\n    </div>\n    <div class=\"col-sm-2\">\n      Diseases <span class=\"badge badge-pill badge-success\">{{diseaseCount}}</span>\n    </div>\n    <div class=\"col-sm-2\">\n      Genes <span class=\"badge badge-pill badge-success\">{{geneCount}}</span>\n    </div>\n  </div>\n  <div class=\"row pb-1\">\n\n    <div class=\"col-sm-3\">\n        <a class=\"btn btn-primary btn-block px-0 mb-1\" style=\"white-space: normal;\" [routerLink]=\"['/pheno-marks',{dataToAdd:'a'}]\" routerLinkActive=\"active\">\n        Pheno Marks</a>\n        <p class=\"small text-center\">paste epicrysis and get phenotypes</p>\n    </div>\n    <div class=\"col-sm-3\">\n        <a class=\"btn btn-primary btn-block px-0 mb-1\" style=\"white-space: normal;\" routerLink=\"/similar-diseases\" routerLinkActive=\"active\">\n        Disease Input</a>\n        <p class=\"small text-center\">input disease and get phenotypes</p>\n    </div>\n    <div class=\"col-sm-3\">\n        <a class=\"btn btn-primary btn-block px-0 mb-1\" style=\"white-space: normal;\" routerLink=\"/pheno-panel\" routerLinkActive=\"active\">\n        Pheno Panel</a>\n        <p class=\"small text-center\">curate or enter phenotypes</p>\n    </div>\n    <!--<div class=\"col-6 col-sm-2\">\n        <a class=\"btn btn-primary btn-block px-0 mb-1\" style=\"white-space: normal;\" routerLink=\"/pheno-panel\" routerLinkActive=\"active\">Genes</a>\n    </div>\n    <div class=\"col-6 col-sm-2\">\n        <a class=\"btn btn-primary btn-block px-0 mb-1\" style=\"white-space: normal;\" routerLink=\"/pheno-panel\" routerLinkActive=\"active\">Diseases</a>\n    </div>-->\n  </div>\n  <router-outlet></router-outlet>\n</div>\n  ",
        }), 
        __metadata('design:paramtypes', [phenotype_service_1.PhenotypeService])
    ], AppComponent);
    return AppComponent;
}());
exports.AppComponent = AppComponent;
//# sourceMappingURL=app.component.js.map