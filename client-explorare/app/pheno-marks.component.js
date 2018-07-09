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
require('rxjs/operator/map');
require('rxjs/Rx');
var phenotype_service_1 = require('./phenotype.service');
var PhenoMarksComponent = (function () {
    function PhenoMarksComponent(phenotypeService) {
        var _this = this;
        this.phenotypeService = phenotypeService;
        this.EXAMPLE_MICRO = "Warburg et al. (1993) used the designation Micro syndrome for an autosomal recessive syndrome comprising microcephaly, microcornea, congenital cataract, mental retardation, optic atrophy, and hypogenitalism. They described an affected brother and sister and their male cousin. The sibs were offspring of a consanguineous Pakistani marriage; the parents of the cousin denied consanguinity. Agenesis of the corpus callosum, prominent root of the nose, large anteverted ears, facial hypertrichosis, small pupils with posterior synechiae, hypotonia, mild to moderate spastic palsy with hip dislocations, and hormonal dysfunction, presumably of hypothalamic origin, were other features. The children were almost blind, whether or not the cataracts had been operated on. The electroretinographic responses indicated dysfunction of both retinal rods and cones, and the visual evoked potentials confirmed optic nerve atrophy. The children were late walkers and were incontinent of urine and stools. In the differential diagnosis, Warburg et al. (1993) considered COFS syndrome (214150), CAMAK/CAMFAK syndromes (212540), Martsolf syndrome (212720), lethal Rutledge syndrome (270400), and lethal Neu-Laxova syndrome (256520).";
        this.query = "";
        this.output = "";
        phenotypeService.parsedText.subscribe(function (value) { return _this.output = value; });
    }
    PhenoMarksComponent.prototype.example = function () {
        this.query = this.EXAMPLE_MICRO;
    };
    PhenoMarksComponent.prototype.postItems = function () {
        this.output = "\n      <div class=\"spinner\">\n        <div class=\"bounce1\"></div>\n        <div class=\"bounce2\"></div>\n        <div class=\"bounce3\"></div>\n      </div>";
        this.phenotypeService.parseText(this.query);
    };
    PhenoMarksComponent = __decorate([
        core_1.Component({
            selector: 'my-pheno-marks',
            template: "\n  <div class=\"row pb-1\">\n    <div class=\"col-sm-9\">\n      <textarea [(ngModel)]=\"query\" rows=\"5\" class=\"border-no width-full\" placeholder=\"Paste epicrysis here\" (change)=postItems()></textarea>\n    </div>\n    <div class=\"col-sm-3\">\n      <a class=\"btn btn-info width-full mb-1\" style=\"color:white;\" (click)=example()>example</a>\n      <a class=\"btn btn-success width-full mb-1\" style=\"color:white;\" (click)=postItems()>submit</a>\n    </div>\n  </div>\n  <div class=\"row pb-1\">\n    <div [innerHTML]=\"output\" style=\"clear: both;\"></div>\n  </div>\n\n  ",
        }), 
        __metadata('design:paramtypes', [phenotype_service_1.PhenotypeService])
    ], PhenoMarksComponent);
    return PhenoMarksComponent;
}());
exports.PhenoMarksComponent = PhenoMarksComponent;
//# sourceMappingURL=pheno-marks.component.js.map