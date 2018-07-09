import { OpaqueToken } from "@angular/core";

export let APP_CONFIG = new OpaqueToken("app.config");

export interface IAppConfig {
    apiEndpoint: string;
}

export const AppConfig: IAppConfig = {
    //apiEndpoint: "http://localhost:8080"
    apiEndpoint: "http://178.216.201.71:8080"
};
