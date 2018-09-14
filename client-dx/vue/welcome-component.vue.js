const Welcome = {
  template:`
    <v-app>
    <v-layout row justify-center>
      <v-dialog v-model="waiting_visibility" persistent max-width="500">
        <v-card>
          <!--<v-card-title class="headline">genomicapt uses cookies</v-card-title>-->
          <v-card-text>signing in</v-card-text>
          <v-card-actions>
            <v-layout column>
              <v-layout row>
                <v-progress-circular
                    :width="1"
                    :size="70"
                    color="grey"
                    indeterminate
                  ></v-progress-circular>
              </v-layout>
            <v-layout>
          </v-card-actions>
        </v-card>
      </v-dialog>
    </v-layout>
    </v-app>`,
    data: function () {
      return {
        waiting_visibility: false
      }
    },
    methods: {
            go: function () {
                this.$router.push("/login")
            }
        }
}
