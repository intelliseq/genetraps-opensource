const sampleModule = {
	namespaced: true,
	state: {
		sample: {properties: {"key1": "value1","key2": "value2"}}
  },
	mutations: {
	},
  actions: {
		createSample({commit}) {
			logger.debug("vue.vuex.sample.createSample")
			request({
				waitingText: "Creating sample",
				endpoint: "sample/new",
				callback: function(data){console.log(data)}
			})
		},
		getSamples({commit}) {
			logger.debug("vue.vuex.sample.getSamples")
			request({
				waitingText: "Fetching samples information",
				endpoint: "user/privileges",
				callback: function(data){console.log(data)}
			})
		},
		showSamples({commit}) {
			logger.debug("vue.vuex.sample.showSamples")
			store.dispatch("sample/getSamples")
		},
	}
}
