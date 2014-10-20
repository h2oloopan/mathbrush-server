hello = require '../addons/build/Release/hello'

exports.bind = (app) ->
	app.get '/api/hello', (req, res) ->
		output = hello.hello()
		res.send 200, output