const loginComponent = {
  template: `
      <v-layout>
          <v-flex>
            <v-dialog v-model="loginDialogVisibility" persistent max-width="500px">
            <v-card class="elevation-6">
              <v-toolbar dark color="primary">
                <v-toolbar-title>Sign in</v-toolbar-title>
              </v-toolbar>
              <v-card-text>
                <v-form>
                  <v-text-field v-model="login" prepend-icon="person" name="login" label="Login" type="text"></v-text-field>
                  <v-text-field v-model="password" prepend-icon="lock" name="password" label="Password" id="password" type="password"></v-text-field>
                </v-form>
              </v-card-text>
              <v-card-actions>
                <v-spacer></v-spacer>
                <v-btn v-on:click.prevent="getToken" color="primary">Login</v-btn>
              </v-card-actions>
            </v-card>
            </v-dialog>
          </v-flex>
          </v-layout>
  `,
  methods: {
          getToken: function (event) {
              logger.debug("vue.login.getToken")
              this.loginDialogVisibility = false
              var credentials = {login: this.login, password: this.password}
              store.dispatch('security/loginWithCredentials', credentials)
          }
      },
      data: function () {
        return {
          login: '',
          password: '',
          loginDialogVisibility: false
        }
      },
      mounted() {
        logger.debug("vue.login.mounted")
        this.loginDialogVisibility = true
      }
}
