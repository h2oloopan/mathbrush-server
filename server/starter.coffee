path = require 'path'
fs = require 'fs'
folder = path.resolve 'apis'

exports.start = (app) ->
	#bind apis
	files = fs.readdirSync folder
	(require path.join(folder, api)).bind app for api in files when path.extname(api) == '.js'