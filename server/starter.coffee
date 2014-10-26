path = require 'path'
fs = require 'fs'
folder = path.resolve 'apis'
index = 'pages/index.html'

exports.start = (app) ->
	#bind apis
	files = fs.readdirSync folder
	(require path.join(folder, api)).bind app for api in files when path.extname(api) == '.js'

	#bind the only index page we need, everything else will be done front-endly
	app.get '/', (req, res) ->
		res.sendfile index, 
			root: __dirname