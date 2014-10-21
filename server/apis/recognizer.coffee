recognizer = require '../addons/build/Release/recognizer'


logger =
	log: (msg) ->
		console.log msg

recognizer.initialize logger

exports.bind = (app) ->
	app.get '/api/recognize/recognize', (req, res) ->
		input = [
			[
				[1, 2]
				[2, 4]
				[3, 6]
				[4, 8]
				[5, 10]
			]
			[
				[7, 10]
				[9, 13]
				[15, 63]
			]
		]


		result = recognizer.recognize 'pan', input
		res.send 200, result