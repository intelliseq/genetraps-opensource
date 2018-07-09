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
var app_config_1 = require('./app.config');
var PhenotypeService = (function () {
    function PhenotypeService(http, config) {
        var _this = this;
        this.http = http;
        this.config = config;
        /*** parsed text ***/
        this.parsedText = new BehaviorSubject_1.BehaviorSubject("");
        /*** phenotypes ***/
        this.phenotypes = new BehaviorSubject_1.BehaviorSubject([]);
        this.diseases = new BehaviorSubject_1.BehaviorSubject([]);
        this.genes = new BehaviorSubject_1.BehaviorSubject([]);
        this.phenotypes.subscribe(function (value) { return _this.getDiseasesByPhenotypes(value); });
        this.phenotypes.subscribe(function (value) { return _this.getGenesByPhenotypes(value); });
    }
    /*** getters & setters ***/
    PhenotypeService.prototype.addPhenotype = function (phenotype) {
        var updatedPhenotypes = this.phenotypes.value;
        updatedPhenotypes.push(phenotype);
        this.phenotypes.next(updatedPhenotypes);
        //this.phenotypes;
    };
    /*** requests ***/
    PhenotypeService.prototype.parseText = function (query) {
        var _this = this;
        var data = JSON.stringify({ "query": query });
        var headers = new http_1.Headers();
        headers.append('Content-Type', 'application/json');
        //console.log(this.query);
        this.http.post(this.config.apiEndpoint + '/parse-text', data, { headers: headers })
            .map(function (response) { return response.json(); })
            .subscribe(function (result) {
            _this.parsedText.next(result.markedText);
            _this.phenotypes.next(result.hpoTerms);
            console.log(result.hpoTerms);
        });
    };
    PhenotypeService.prototype.getPhenotypesByDisease = function (disease) {
        var _this = this;
        var parsedDisease = disease.replace(/ /g, '%20');
        var url = this.config.apiEndpoint + "/get-phenotypes-by-disease?disease=" + parsedDisease;
        this.http.get(url)
            .map(function (response) { return response.json(); })
            .subscribe(function (result) {
            _this.phenotypes.next(result.results[0]);
            console.log(result.hpoTerms);
        });
    };
    PhenotypeService.prototype.getDiseasesByPhenotypes = function (phenotypes) {
        var _this = this;
        //console.log("disease request");
        var data = "{ \"hpoTerms\": " + JSON.stringify(phenotypes) + " }";
        var headers = new http_1.Headers();
        headers.append('Content-Type', 'application/json');
        //console.log(this.query);
        this.http.post(this.config.apiEndpoint + '/get-diseases', data, { headers: headers })
            .map(function (response) { return response.json(); })
            .subscribe(function (result) {
            _this.diseases.next(result);
            //console.log(result);
            console.log("disease request completed");
        });
    };
    PhenotypeService.prototype.getGenesByPhenotypes = function (phenotypes) {
        var _this = this;
        //console.log("disease request");
        var data = "{ \"hpoTerms\": " + JSON.stringify(phenotypes) + " }";
        var headers = new http_1.Headers();
        headers.append('Content-Type', 'application/json');
        //console.log(this.query);
        this.http.post(this.config.apiEndpoint + '/get-genes', data, { headers: headers })
            .map(function (response) { return response.json(); })
            .subscribe(function (result) {
            _this.genes.next(result);
            //console.log(result);
            console.log("gene request completed");
        });
    };
    PhenotypeService = __decorate([
        core_1.Injectable(),
        __param(1, core_1.Inject(app_config_1.APP_CONFIG)), 
        __metadata('design:paramtypes', [http_1.Http, Object])
    ], PhenotypeService);
    return PhenotypeService;
}());
exports.PhenotypeService = PhenotypeService;
//# sourceMappingURL=phenotype.service.js.map