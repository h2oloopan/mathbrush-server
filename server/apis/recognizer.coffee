recognizer = require '../addons/build/Release/recognizer'


logger =
	log: (msg) ->
		console.log msg

recognizer.initialize logger

exports.bind = (app) ->
	app.post '/api/recognize/recognize', (req, res) ->
		result = recognizer.recognize 'pan', req.body
		res.send 200, result