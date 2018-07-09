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
var __param = (this && this.__param) || function (paramIndex, decorator) {
    return function (target, key) { decorator(target, key, paramIndex); }
};
var core_1 = require('@angular/core');
require('rxjs/operator/map');
require('rxjs/Rx');
var phenotype_service_1 = require('./phenotype.service');
var app_config_1 = require('./app.config');
var SimilarDiseasesComponent = (function () {
    function SimilarDiseasesComponent(phenotypeService, config) {
        this.phenotypeService = phenotypeService;
        this.config = config;
        this.diseaseAutocompleteAdress = this.config.apiEndpoint + "/disease-autocomplete?firstLetters=:keyword&resultsCount=20";
    }
    SimilarDiseasesComponent.prototype.getPhenotypes = function (disease) {
        this.phenotypeService.getPhenotypesByDisease(disease);
    };
    SimilarDiseasesComponent = __decorate([
        core_1.Component({
            selector: 'my-similar-diseases',
            template: "\n  <div class=\"row pb-1\">\n    <div class=\"col-xs-12\">\n    <h4>Input disease</h4>\n        <input  ng2-auto-complete\n          [(ngModel)]=\"diseaseName\"\n          (ngModelChange)=\"getPhenotypes($event)\"\n          placeholder=\"Enter Disease Name\"\n          [source]=\"diseaseAutocompleteAdress\"\n          path-to-data=\"results\"\n          min-chars=\"0\" />\n    </div>\n  </div>\n  <div class=\"row pb-1\">\n    <div class=\"col-xs-12\">\n    <h4>{{diseaseName}}</h4>\n    </div>\n  </div>\n  ",
            styles: ["\n    ng2-auto-complete, input {\n      display: block; border: 0px solid #ccc; width: 600px;\n    }\n  "],
        }),
        __param(1, core_1.Inject(app_config_1.APP_CONFIG)), 
        __metadata('design:paramtypes', [phenotype_service_1.PhenotypeService, Object])
    ], SimilarDiseasesComponent);
    return SimilarDiseasesComponent;
}());
exports.SimilarDiseasesComponent = SimilarDiseasesComponent;
//# sourceMappingURL=similar-diseases.component.js.map