var loadGoogleAnalytics = function(){
  logger("DEBUG", "loadGoogleAnalytics");

  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-125435882-1', 'auto');
  ga('send', 'pageview');

    ga('set', 'page', router.currentRoute.path);
    ga('send', 'pageview');

    router.afterEach(( to, from ) => {
      ga('set', 'page', to.path);
      ga('send', 'pageview');
      logger("DEBUG", "vue.app.router.afterEach " + to.path);
    });

 }