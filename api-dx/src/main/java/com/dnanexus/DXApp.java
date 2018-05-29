// Copyright (C) 2013-2016 DNAnexus, Inc.
//
// This file is part of dx-toolkit (DNAnexus platform client libraries).
//
//   Licensed under the Apache License, Version 2.0 (the "License"); you may
//   not use this file except in compliance with the License. You may obtain a
//   copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
//   WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
//   License for the specific language governing permissions and limitations
//   under the License.

package com.dnanexus;

import com.dnanexus.DXHTTPRequest.RetryStrategy;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonInclude.Include;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.JsonNode;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

import java.util.List;
import java.util.Map;

/**
 * An applet (an executable data object).
 *
 * <p>
 * Although these bindings will allow you to create a simple applet from scratch, we encourage you
 * to use the command-line tool <code>dx build</code> instead. See the <a
 * href="https://wiki.dnanexus.com/API-Specification-v1.0.0/Applets-and-Entry-Points">API
 * documentation for applets</a>.
 * </p>
 *
 * <p>
 * The {@link #describe()} method returns a {@link RunSpecification} object that does not have the
 * "code" field populated; that is, {@link RunSpecification#getCode()} will return {@code null}.
 * </p>
 */
public class DXApp extends DXApplet implements DXExecutable<DXJob> {

    /**
     * Contains metadata for an applet.
     */
    public static class Describe extends DXDataObject.Describe {
        @JsonProperty
        private String title;
        @JsonProperty
        private String summary;
        @JsonProperty
        private String description;
        @JsonProperty
        private List<InputParameter> inputSpec;
        @JsonProperty
        private List<OutputParameter> outputSpec;
        @JsonProperty
        private RunSpecification runSpec;
        @JsonProperty
        private String dxapi;

        private Describe() {
            super();
        }

        /**
         * Returns the applet description.
         *
         * @return applet description
         */
        public String getDescription() {
            Preconditions
                    .checkState(this.description != null,
                            "description is not available because it was not retrieved with the describe call");
            return description;
        }

        /**
         * Returns the API version that the applet code is to run under.
         *
         * @return API version
         */
        public String getDXAPIVersion() {
            Preconditions
                    .checkState(this.dxapi != null,
                            "dxapi version is not available because it was not retrieved with the describe call");
            return dxapi;
        }

        /**
         * Returns the applet's input specification.
         *
         * @return applet input specification
         */
        public List<InputParameter> getInputSpecification() {
            Preconditions
                    .checkState(this.inputSpec != null,
                            "input specification is not available because it was not retrieved with the describe call");
            return ImmutableList.copyOf(this.inputSpec);
        }

        /**
         * Returns the applet's output specification.
         *
         * @return applet output specification
         */
        public List<OutputParameter> getOutputSpecification() {
            Preconditions
                    .checkState(this.outputSpec != null,
                            "output specification is not available because it was not retrieved with the describe call");
            return ImmutableList.copyOf(this.outputSpec);
        }

        /**
         * Returns the applet's run specification.
         *
         * @return applet run specification
         */
        public RunSpecification getRunSpecification() {
            Preconditions
                    .checkState(this.runSpec != null,
                            "run specification is not available because it was not retrieved with the describe call");
            return runSpec;
        }

        /**
         * Returns the applet summary.
         *
         * @return applet summary
         */
        public String getSummary() {
            Preconditions.checkState(this.summary != null,
                    "summary is not available because it was not retrieved with the describe call");
            return summary;
        }

        /**
         * Returns the applet title.
         *
         * @return applet title
         */
        public String getTitle() {
            Preconditions.checkState(this.title != null,
                    "title is not available because it was not retrieved with the describe call");
            return title;
        }
    }

    /**
     * Deserializes a DXApplet from JSON containing a DNAnexus link.
     *
     * @param value JSON object map
     *
     * @return data object
     */
    @JsonCreator
    private static DXApp create(Map<String, Object> value) {
        checkDXLinkFormat(value);
        // TODO: how to set the environment?
        return DXApp.getInstance((String) value.get("$dnanexus_link"));
    }

    /**
     * Returns a {@code DXApplet} associated with an existing applet.
     *
     * @throws NullPointerException If {@code appletId} is null
     */
    public static DXApp getInstance(String appletId) {
        return new DXApp(appletId, null);
    }

    /**
     * Returns a {@code DXApplet} associated with an existing applet in a particular project or
     * container.
     *
     * @throws NullPointerException If {@code appletId} or {@code container} is null
     */
    public static DXApp getInstance(String appletId, DXContainer project) {
        return new DXApp(appletId, project, null, null);
    }

    /**
     * Returns a {@code DXApplet} associated with an existing applet in a particular project using
     * the specified environment, with the specified cached describe output.
     *
     * <p>
     * This method is for use exclusively by bindings to the "find" routes when describe hashes are
     * returned with the find output.
     * </p>
     *
     * @throws NullPointerException If any argument is null
     */
    static DXApp getInstanceWithCachedDescribe(String appletId, DXContainer project,
                                               DXEnvironment env, JsonNode describe) {
        return new DXApp(appletId, project, Preconditions.checkNotNull(env,
                "env may not be null"), Preconditions.checkNotNull(describe,
                "describe may not be null"));
    }

    /**
     * Returns a {@code DXApplet} associated with an existing applet in a particular project using
     * the specified environment.
     *
     * @throws NullPointerException If {@code appletId} or {@code container} is null
     */
    public static DXApp getInstanceWithEnvironment(String appletId, DXContainer project,
                                                   DXEnvironment env) {
        return new DXApp(appletId, project, Preconditions.checkNotNull(env,
                "env may not be null"), null);
    }

    /**
     * Returns a {@code DXApplet} associated with an existing applet using the specified
     * environment.
     *
     * @throws NullPointerException If {@code appletId} is null
     */
    public static DXApp getInstanceWithEnvironment(String appletId, DXEnvironment env) {
        return new DXApp(appletId, Preconditions.checkNotNull(env, "env may not be null"));
    }

    private DXApp(String appletId, DXContainer project, DXEnvironment env, JsonNode describe) {
        super(appletId, "app", project, env, describe);
    }

    private DXApp(String appletId, DXEnvironment env) {
        super(appletId, "app", env);
    }


}
