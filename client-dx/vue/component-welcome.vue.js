const welcomeComponent = {
  template:`
  <v-container fluid grid-list-sm>
  <v-layout row wrap pa-1>
    <v-flex xs12 sm6 md4 pa-2>
    <v-card color="cyan darken-1" class="white--text">
      <v-layout align-center>
        <v-flex xs4 pa-2>
          <v-img src="images/abstract6.png" height="50%""></v-img>
        </v-flex>
        <v-flex xs8>
          <v-card-title primary-title>
          <div>
            <div class="headline">Get samples</div>
            <div>Create sample, describe sample, upload fastq</div>
            <v-btn v-on:click.prevent="getSamples" dark color="grey darken-3">Create now</v-btn>
          </div>
          </v-card-title>
        </v-flex>
      </v-layout>
      </v-card>
    </v-flex>
    <v-flex xs12 sm6 md4 pa-2>
    <v-card color="cyan darken-3" class="white--text">
      <v-layout align-center>
        <v-flex xs4 pa-2>
          <v-img src="images/abstract3.png" height="50%""></v-img>
        </v-flex>
        <v-flex xs8>
          <v-card-title primary-title>
          <div>
            <div class="headline">Create sample</div>
            <div>Create sample, describe sample, upload fastq</div>
            <v-btn v-on:click.prevent="createSample" dark color="grey darken-3">Create now</v-btn>
          </div>
          </v-card-title>
        </v-flex>
      </v-layout>
      </v-card>
    </v-flex>
    <v-flex xs12 sm6 md4 pa-2>
    <v-card color="cyan darken-2" class="white--text">
      <v-layout align-center>
        <v-flex xs4 pa-2>
          <v-img src="images/abstract4.png" height="50%""></v-img>
        </v-flex>
        <v-flex xs8>
          <v-card-title primary-title>
          <div>
            <div class="headline">Create sample</div>
            <div>Create sample, describe sample, upload fastq</div>
            <v-btn v-on:click.prevent="createSample" dark color="grey darken-3">Create now</v-btn>
          </div>
          </v-card-title>
        </v-flex>
      </v-layout>
      </v-card>
    </v-flex>
    <v-flex xs12 sm6 md4 pa-2>
    <v-card color="cyan darken-1" class="white--text">
      <v-layout align-center>
        <v-flex xs4 pa-2>
          <v-img src="images/abstract5.png" height="50%""></v-img>
        </v-flex>
        <v-flex xs8>
          <v-card-title primary-title>
          <div>
            <div class="headline">Create sample</div>
            <div>Create sample, describe sample, upload fasta</div>
            <v-btn v-on:click.prevent="createSample" dark color="grey darken-3">Create now</v-btn>
          </div>
          </v-card-title>
        </v-flex>
      </v-layout>
      </v-card>
    </v-flex>
  </v-layout>
  </v-container>
  <!--<v-layout row wrap pa-4>
        <v-flex xs12 sm6 md3>
          <v-card color="blue-grey darken-2" class="white--text">
            <v-card-title primary-title>
              <div class="headline">Create sample</div>
            </v-card-title>
            <v-card-actions>
              <v-btn v-on:click.prevent="createSample" flat dark>Create now</v-btn>
            </v-card-actions>
          </v-card>
        </v-flex>
    </v-layout>-->

    `,
    data: function () {
      return {

      }
    },
    methods: {
      createSample: function (event) {
        logger.debug("vue.welcome.createSample")
        store.dispatch('sample/createSample')
      },
      getSamples: function (event) {
        logger.debug("vue.welcome.getSamples")
        store.dispatch('user/getSamples')
      },
    },
    created: function () {
      logger.debug("vue.welcome.created")
    }
}
