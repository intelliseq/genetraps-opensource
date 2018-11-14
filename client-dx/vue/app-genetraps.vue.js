const router = new VueRouter({
  mode: 'history',
  routes: [
    { path: '/', component: welcomeComponent },
    { path: '/login', component: loginComponent },
    { path: '/cookies', component: cookiesComponent }
    { path: '/samples', component: samplesComponent }
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
    logger.debug("vue.app.created")
  },
  mounted() {

    /* first we check for cookies */
    /* then we check for google anlytics */
    /* finally we check for refresh token */
    logger.debug("vue.app.mounted")
    if(this.$cookies.get("cookies_necessary") == null) {
      logger.debug("vue.app no agreement for cookies")
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
    'wait-component': waitComponent,
    'toolbar-component': toolbarComponent,
    'footer-component': footerComponent,
    'samples-component': samplesComponent,
    'left-toolbar-component': leftToolbarComponent
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
