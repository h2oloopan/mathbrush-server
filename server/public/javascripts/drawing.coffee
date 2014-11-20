url_recognize = 'api/recognize/recognize' 


class Drawing
	constructor: (canvas, dots, preview, mathML, latex) ->
		@canvas = canvas[0]
		@dots = dots[0]
		@previewHolder = preview
		@mathMLHolder = mathML
		@latexHolder = latex
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
		thiz = @
		$(@dots).hide()
		$(@canvas).show()
		data = JSON.stringify @strokes
		$.ajax
			type: 'POST'
			contentType: 'application/json; charset=utf-8'
			url: url_recognize
			data: data
			dataType: 'json'
		.done (result) ->
			mathML = result[0]
			latex = result[1]
			thiz.display mathML, latex
		.fail (response) ->
			alert response.responseText
	display: (mathML, latex) ->
		@mathMLHolder.text mathML
		@latexHolder.text latex
		#display true math
		@previewHolder.html '$' + latex + '$'
		MathJax.Hub.Queue ['Typeset', MathJax.Hub, @previewHolder[0]]
	drawDots: (button) ->
		$(@canvas).hide()
		$(@dots).show()
		ctx = @dots.getContext '2d'
		ctx.fillStyle = '#EEEEEE'
		ctx.fillRect 0, 0, @dots.width, @dots.height

		#draw the strokes
		getRandomColor = ->
			letters = '0123456789ABCDEF'.split ''
			color = '#'
			for i in [0...6]
				color += letters[Math.floor(Math.random() * 16)]
			return color

		for stroke in @strokes
			color = getRandomColor()
			for point in stroke
				x = point[0]
				y = point[1]
				#console.log x + ' ' + y
				ctx.fillStyle = color
				ctx.fillRect x, y, 5, 5

	clean: ->
		@clicked = 0
		@stroke = []
		@strokes = []
		@ctx.fillStyle = '#FFFFFF'
		@ctx.fillRect 0, 0, @canvas.width, @canvas.height
		#clean ui
		@previewHolder.html ''
		@latexHolder.html ''
		@mathMLHolder.html ''


