recognizer = require '../addons/build/Release/recognizer'

log = (msg) ->
	console.log msg

recognizer.initialize log

exports.bind = (app) ->
	app.get '/api/recognize/recognize', (req, res) ->
		result = recognizer.recognize 'pan', [1, 2, 3, 4]
		res.send 200, result