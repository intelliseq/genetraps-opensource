QUnit.test( "hello test", function( assert ) {
  assert.ok( 1 == "1", "Passed!" );
});

QUnit.test('"Hello Vue.js!" is displayed', function (assert) {
  assert.equal(vm.$data.message, 'Hello Vue.js!', "Passed!");
  //assert.equal(document.querySelector('h1').textContent, 'Hello Vue.js!');
});
