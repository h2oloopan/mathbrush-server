class Drawing
	constructor: (canvas) ->
		@canvas = canvas[0]
	draw: ->
		#this is the initialize function
		ctx = @canvas.getContext '2d'
		ctx.fillStyle = '#FFFFFF'
		ctx.fillRect 0, 0, @canvas.width, @canvas.height