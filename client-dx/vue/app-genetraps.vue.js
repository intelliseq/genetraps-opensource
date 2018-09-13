/* event hub */
Vue.prototype.$hub = new Vue();

const router = new VueRouter({
  mode: 'history',
  routes: [
    { path: '/', component: Welcome },
    { path: '/login', component: Login },
    { path: '/cookies', component: CookiesComponent }
  ]
})

new Vue({
	router,
  el: '#genetraps-app',
  data: {
    access_token: ""
  },
  created() {
    console.log("LOG: Vue.App.created()")
    this.$hub.$on('login', (response) => {
      this.$cookies.set("refresh_token", response.refresh_token, "7D")
      this.access_token = response.access_token
      console.log("LOG: Vue.App.hub.login " + this.access_token)

    });
  },
  mounted() {
    console.log("LOG: Vue.App.mounted()")
    if(this.$cookies.get("cookies_necessary") == null) {
      console.log("LOG: Vue.App.mounted() no agreement for cookies")
      this.$router.push("/cookies")
    } else {
      if(this.$cookies.get("cookies_statistics") != null) {
        loadGoogleAnalytics()
      }
    }

    //if(this.$cookies.get("refresh_token") == null) {
    //  this.$router.push("/login")
    //} else {
    //  this.loginWithRefreshToken()
    //}
  },
  methods: {
    loginWithRefreshToken() {
      var reqData = {
        "grant_type": "refresh_token",
        "client_id": "web_app",
        "refresh_token": this.$cookies.get("refresh_token")
      }
      axios({
        method: 'post', //you can set what request you want to be
        url: 'http://genetraps.intelliseq.pl:8088/oauth/token',
        withCredentials: true,
        crossdomain: true,
        data: Object.keys(reqData).map(function(key) {
          return encodeURIComponent(key) + '=' + encodeURIComponent(reqData[key])
        }).join('&'),
        headers: {
          'Authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
          'Content-Type': 'application/x-www-form-urlencoded' },
      })
      .then(response => {
        console.log("LOG: Vue.App.loginWithRefreshToken() succesfull login")
        this.$hub.$emit('login', response.data);
        //console.log(response.data.access_token)
      })
      .catch(e => {
        console.log(e)
        this.error = true
        this.dialog_visibility = false
        this.error_message = "Error: " + e.response.data.error
        this.error_description = "Error description: " + e.response.data.error_description
        this.alert_visibility = true
      })
    }

  },
  components: {
    'login-component': Login
  }
})

// google analytics and router cooperation
//ga('set', 'page', router.currentRoute.path);
//ga('send', 'pageview');

/*router.afterEach(( to, from ) => {
  ga('set', 'page', to.path);
  ga('send', 'pageview');
  console.console.log("Vue.App.router.afterEach() " + to.path);
});*/
