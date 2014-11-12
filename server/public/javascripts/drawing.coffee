class Drawing
	constructor: (canvas) ->
		@canvas = canvas[0]
	draw: ->
		#this is the initialize function
		ctx = @canvas.getContext '2d'
		ctx.fillStyle = '#FFFFFF'
		ctx.fillRect 0, 0, @canvas.width, @canvas.height
		ctx.lineWidth = 3
		ctx.lineCap = 'round'
		#bind
		getMousePos = (canvas, e) ->
			rect = canvas.getBoundingClientRect()
			return pos =
				x: e.clientX - rect.left
				y: e.clientY - rect.top

		drawMouse = (canvas) ->
			clicked = 0
			start = (e) ->
				clicked = 1
				ctx.beginPath()
				pos = getMousePos canvas, e
				ctx.moveTo pos.x, pos.y

			move = (e) ->
				if clicked
					pos = getMousePos canvas, e
					ctx.lineTo pos.x, pos.y
					ctx.stroke()

			stop = (e) ->
				clicked = 0

			$(canvas).on 'mousedown', start
			$(canvas).on 'mousemove', move
			$(canvas).on 'mouseup', stop

		drawMouse @canvas
	recognize: ->
		alert 'recognize'
	clean: ->
		alert 'clean'