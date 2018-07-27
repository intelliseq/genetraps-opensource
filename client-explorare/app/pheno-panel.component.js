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
var BehaviorSubject_1 = require('rxjs/BehaviorSubject');
var http_1 = require('@angular/http');
require('rxjs/operator/map');
require('rxjs/Rx');
var phenotype_service_1 = require('./phenotype.service');
var phenotype_1 = require('./phenotype');
var app_config_1 = require('./app.config');
var PhenoPanelComponent = (function () {
    function PhenoPanelComponent(phenotypeService, config, http, ref) {
        var _this = this;
        this.phenotypeService = phenotypeService;
        this.config = config;
        this.http = http;
        this.ref = ref;
        this.phenotypeAutocompleteAdress = this.config.apiEndpoint + "/phenotype-autocomplete?firstLetters=:keyword&resultsCount=10";
        /* dialog */
        this.showDialog = false;
        this.showSpinner = false;
        this.dialogPhenotype = new BehaviorSubject_1.BehaviorSubject("");
        this.dialogObject = null;
        this.dialogPhenotypeOriginId = "";
        phenotypeService.phenotypes.subscribe(function (value) {
            console.log("constructorSubscribePhenotypes");
            _this.phenotypes = value;
        });
        phenotypeService.diseases.subscribe(function (value) { return _this.diseases = value; });
        phenotypeService.genes.subscribe(function (value) {
            _this.genes = value;
            _this.getGeneNames(value);
        });
        this.dialogPhenotype.subscribe(function (value) {
            console.log("constructorSubscribeDialogPhenotype: " + value);
            _this.getPhenotype(value);
        });
    }
    PhenoPanelComponent.prototype.dialog = function (id) {
        console.log("dialog: " + id);
        this.dialogPhenotype.next(id);
    };
    PhenoPanelComponent.prototype.close = function () {
        this.showDialog = false;
    };
    PhenoPanelComponent.prototype.getPhenotype = function (id) {
        var _this = this;
        this.showSpinner = true;
        //console.log("get-phenotype-1: " + id);
        var url = this.config.apiEndpoint + "/get-phenotype-by-id?phenotypeid=" + id;
        //console.log("get-phenotype-2: " + id);
        //this.dialogObject = null;
        //console.log("get-phenotype-3: " + id);
        if (id != "") {
            //console.log("get-phenotype-4: " + id);
            this.http.get(url)
                .map(function (response) { return response.json(); })
                .subscribe(function (result) {
                console.log("get-phenotype-subscribe: " + result.phenotype.id);
                _this.showSpinner = false;
                _this.dialogObject = result;
                console.log(result);
                _this.showDialog = true;
            });
        }
        console.log("get-phenotype-5: " + id);
    };
    PhenoPanelComponent.prototype.go = function (id) {
        console.log("go: " + id);
        this.dialogPhenotypeOriginId = id;
    };
    PhenoPanelComponent.prototype.usePhenotype = function (id, name) {
        console.log("use-phenotype: " + id);
        this.addPhenotype(new phenotype_1.Phenotype(id, name));
        this.remove(this.dialogPhenotypeOriginId);
        this.dialogPhenotypeOriginId = "";
        this.showSpinner = false;
        this.showDialog = false;
        this.dialogObject = null;
    };
    PhenoPanelComponent.prototype.addPhenotype = function (phenotype) {
        console.log("add-phenotype: " + phenotype.id);
        if (phenotype.name != null) {
            this.phenotypeService.addPhenotype(phenotype);
        }
    };
    PhenoPanelComponent.prototype.remove = function (id) {
        console.log("remove: " + id);
        for (var i = 0; i < this.phenotypes.length; i++) {
            if (this.phenotypes[i].id === id) {
                this.phenotypes.splice(i, 1);
                this.phenotypeService.phenotypes.next(this.phenotypes);
                break;
            }
        }
    };
    PhenoPanelComponent.prototype.getGeneNames = function (genes) {
        var output = "";
        for (var _i = 0, genes_1 = genes; _i < genes_1.length; _i++) {
            var gene = genes_1[_i];
            //console.log("name:")
            //console.log(gene.name)
            output = output + gene.name + "\n"; // 1, "string", false
        }
        this.genesText = output;
    };
    PhenoPanelComponent = __decorate([
        core_1.Component({
            selector: 'my-pheno-panel',
            template: "\n\n\n  <div class=\"row pb-1\">\n  <div class=\"col-sm-3\">\n      add phenotype<input  ng2-auto-complete\n        [(ngModel)]=\"phenotypeName\"\n        placeholder=\"enter phenotype name\"\n        (ngModelChange)=\"addPhenotype($event)\"\n        [source]=\"phenotypeAutocompleteAdress\"\n        value-property-name=null\n        display-property-name=\"name\"\n        path-to-data=\"results\"\n        min-chars=\"0\"\n        style=\"width: 100%\"\n        class=\"pt-1\"/>\n  </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/head.png\" (click)=\"dialog('HP:0000152')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/eye.png\" (click)=\"dialog('HP:0000478')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/connect.png\" (click)=\"dialog('HP:0003549')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/bones.png\" (click)=\"dialog('HP:0000924')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/brain.png\" (click)=\"dialog('HP:0000707')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/digest.png\" (click)=\"dialog('HP:0025031')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/heart.png\" (click)=\"dialog('HP:0001626')\">\n    </div>\n    <div class=\"col-sm-1\">\n    <img src=\"app/images/muscle.png\" (click)=\"dialog('HP:0003011')\">\n    </div>\n  </div>\n  <div *ngIf=\"phenotypes.length > 0\" class=\"row pb-1 border-between\">\n    <div class=\"col-sm-5 pt-2\">\n      <h5>Phenotypes</h5>\n      <ul>\n        <li *ngFor=\"let phenotype of phenotypes\">\n          <span class=\"badge badge-info\">{{phenotype.id}}</span> {{phenotype.name}}\n          <span class=\"badge badge-pill badge-warning\" (click)=\"go(phenotype.id);dialog(phenotype.id)\">more</span>\n          <span class=\"badge badge-pill badge-danger\" (click)=\"remove(phenotype.id)\">remove</span>\n        </li>\n      </ul>\n    </div>\n    <div class=\"col-sm-5 pt-2\">\n      <h5>Diseases <span class=\"badge badge-pill badge-default\">% match</span></h5>\n      <ul>\n        <li *ngFor=\"let disease of diseases\">\n          {{disease.name}} <span class=\"badge badge-pill badge-default\">{{disease.score.toFixed(1)}}</span>\n        </li>\n      </ul>\n    </div>\n    <div class=\"col-sm-2 pt-2\">\n      <button class=\"btn btn-default\" type=\"button\" ngxClipboard [cbContent]=\"genesText\">copy</button>\n\n      <h5>Genes <span class=\"badge badge-pill badge-default\">% match</span></h5>\n\n      <ul>\n        <li *ngFor=\"let gene of genes\">\n          {{gene.name}} <span class=\"badge badge-pill badge-default\">{{gene.score.toFixed(1)}}</span>\n        </li>\n      </ul>\n    </div>\n  </div>\n\n  <div *ngIf=\"showDialog\" class=\"dialog\">\n      <div *ngIf=\"showSpinner\" class=\"spinner\">\n        <div class=\"bounce1\"></div>\n        <div class=\"bounce2\"></div>\n        <div class=\"bounce3\"></div>\n      </div>\n      <h5> Phenotype </h5>\n      <span class=\"badge badge-info\">{{dialogObject.phenotype.id}}</span> {{dialogObject.phenotype.name}}\n      <span class=\"badge badge-pill badge-primary\" (click)=\"usePhenotype(dialogObject.phenotype.id,dialogObject.phenotype.name)\">add</span>\n      <div  *ngIf=\"dialogObject.parent != null\" >\n      <h5 class=\"mt-4\"> Superclass </h5>\n      <span class=\"badge badge-info\">{{dialogObject.parent.id}}</span> {{dialogObject.parent.name}}\n      <span class=\"badge badge-pill badge-info\" (click)=\"dialog(dialogObject.parent.id)\">more</span>\n      <span class=\"badge badge-pill badge-primary\" (click)=\"usePhenotype(dialogObject.parent.id,dialogObject.parent.name)\">add</span>\n      </div>\n      <h5 *ngIf=\"dialogObject.children.length > 0\" class=\"mt-4\"> Subclasses </h5>\n      <ul>\n        <li *ngFor=\"let phenotype of dialogObject.children\">\n          <span class=\"badge badge-info\">{{phenotype.id}}</span> {{phenotype.name}}\n          <span class=\"badge badge-pill badge-info\" (click)=\"dialog(phenotype.id)\">more</span>\n          <span class=\"badge badge-pill badge-primary\" (click)=\"usePhenotype(phenotype.id, phenotype.name)\">add</span>\n        </li>\n      </ul>\n    <button (click)=\"close()\" class=\"btn\">Close</button>\n  </div>\n  <div *ngIf=\"showDialog\" class=\"overlay\" (click)=\"close()\"></div>",
            styles: ["\n  ng2-auto-complete, input {\n    display: block; border: 0px solid #ccc; width: 300px;\n  }\n .border-between > [class*='col-']:before {\n background: #999999;\n bottom: 0;\n content: \" \";\n left: 0;\n position: absolute;\n width: 1px;\n top: 0;\n}\n.border-between > [class*='col-']:first-child:before {\n display: none;\n}\n.overlay {\n  position: fixed;\n  top: 0;\n  bottom: 0;\n  left: 0;\n  right: 0;\n  background-color: rgba(0, 0, 0, 0.5);\n  z-index: 999;\n}\n\n.dialog {\n  z-index: 1000;\n  position: fixed;\n  right: 0;\n  left: 0;\n  top: 20px;\n  margin-right: auto;\n  margin-left: auto;\n  min-height: 200px;\n  width: 90%;\n  max-width: 820px;\n  max-height: calc(100vh - 50px);\n  overflow-y: auto;\n  background-color: #fff;\n  padding: 12px;\n  box-shadow: 0 7px 8px -4px rgba(0, 0, 0, 0.2), 0 13px 19px 2px rgba(0, 0, 0, 0.14), 0 5px 24px 4px rgba(0, 0, 0, 0.12);\n}\n\n@media (min-width: 768px) {\n  .dialog {\n    top: 40px;\n  }\n}\n\n.dialog__close-btn {\n  border: 0;\n  background: none;\n  color: #2d2d2d;\n  position: absolute;\n  top: 8px;\n  right: 8px;\n  font-size: 1.2em;\n}\n"],
        }),
        __param(1, core_1.Inject(app_config_1.APP_CONFIG)), 
        __metadata('design:paramtypes', [phenotype_service_1.PhenotypeService, Object, http_1.Http, core_1.ChangeDetectorRef])
    ], PhenoPanelComponent);
    return PhenoPanelComponent;
}());
exports.PhenoPanelComponent = PhenoPanelComponent;
//# sourceMappingURL=pheno-panel.component.js.map