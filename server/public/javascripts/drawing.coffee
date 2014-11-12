url_recognize = 'api/recognize' 


class Drawing
	constructor: (canvas) ->
		@canvas = canvas[0]
		@ctx = @canvas.getContext '2d'
		@strokes = []
		@stroke = []
		@clicked = 0
	getMousePos: (e) ->
		rect = @canvas.getBoundingClientRect()
		return pos = 
			x: e.clientX - rect.left
			y: e.clientY - rect.top
	start: (e) ->
		p = @getMousePos e
		@clicked = 1
		@ctx.beginPath()
		@ctx.moveTo p.x, p.y
		@stroke = [[p.x, p.y]]
	move: (e) ->
		if @clicked
			p = @getMousePos e
			@ctx.lineTo p.x, p.y
			@ctx.stroke()
			@stroke.push [p.x, p.y]
	stop: (e) ->
		@clicked = 0
		@strokes.push @stroke
	init: ->
		#this is the initialize function
		@ctx.fillStyle = '#FFFFFF'
		@ctx.fillRect 0, 0, @canvas.width, @canvas.height
		@ctx.lineWidth = 3
		@ctx.lineCap = 'round'
		#bind
		$(@canvas).on 'mousedown', (e) =>
			@start e
		$(@canvas).on 'mousemove', (e) =>
			@move e
		$(@canvas).on 'mouseup', (e) =>
			@stop e
	recognize: ->
		data = @strokes
		$.ajax
			type: 'POST'
			contentType: 'application/json; charset=utf-8'
			url: url_recognize
			data: data
			dataType: 'json'
		.done (result) ->
			alert result
		.fail (response) ->
			alert response.responseText
	clean: ->
		@clicked = 0
		@stroke = []
		@strokes = []
		@ctx.fillStyle = '#FFFFFF'
		@ctx.fillRect 0, 0, @canvas.width, @canvas.height