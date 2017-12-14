import Vue from 'vue'
import Hello from '../../src/components/Hello.vue'

describe('Hello.vue', () => {

  it('should render correct contents', () => {
    const vm = new Vue({
      el: document.createElement('div'),
     render: (h) => h(Hello)
    })
    expect(vm.$el.querySelector('h1').textContent).toBe('Welcome to Explorare Client')
  })

  it('has a data function', () => {
    expect(typeof Hello.data).toBe('function')
  })

  it('sets the correct default data', () => {
    const defaultData = Hello.data()
    expect(defaultData.msg).toBe('Welcome to Explorare Client')
  })

  it('should fetch data from api', (done) => {
    const vm = new Vue(Hello).$mount()
    expect(vm.data).toBe('Empty')
    vm.getData()
    setTimeout(function() {
	expect(vm.data.name).toBe('Luke Skywalker'); 
	done()}, 
	2000);
  })

})

