import { OpaqueToken } from "@angular/core";

export let APP_CONFIG = new OpaqueToken("app.config");

export interface IAppConfig {
    apiEndpoint: string;
}

export const AppConfig: IAppConfig = {
    apiEndpoint: "http://localhost:8082"
    //apiEndpoint: "http://genetraps.intelliseq.pl:8082"
};
