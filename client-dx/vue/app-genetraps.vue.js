const router = new VueRouter({
  mode: 'history',
  routes: [
    { path: '/', component: welcomeComponent },
    { path: '/login', component: loginComponent },
    { path: '/cookies', component: cookiesComponent }
  ]
})

const app = new Vue({
	router,
  store,
  el: '#genetraps-app',
  data: {
    access_token: ""
  },
  created() {
    logger("DEBUG", "vue.app.created")
  },
  mounted() {

    /* first we check for cookies */
    /* then we check for google anlytics */
    /* finally we check for refresh token */
    logger("DEBUG", "vue.app.mounted")
    if(this.$cookies.get("cookies_necessary") == null) {
      logger("DEBUG", "vue.app no agreement for cookies")
      this.$router.push("/cookies")
    } else {
      if(this.$cookies.get("cookies_statistics") != null) {
        loadGoogleAnalytics()
      }
      if(this.$cookies.get("refresh_token") == null) {
        this.$router.push("/login")
      } else {
        store.dispatch('security/loginWithRefreshToken')
      }
    }


  },
  components: {
    'wait': waitComponent,
    'toolbar': toolbarComponent
  },
  methods: {

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