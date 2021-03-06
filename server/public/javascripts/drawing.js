// Generated by CoffeeScript 1.7.1
var Drawing, url_recognize;

url_recognize = 'api/recognize/recognize';

Drawing = (function() {
  function Drawing(canvas, dots, preview, mathML, latex) {
    this.canvas = canvas[0];
    this.dots = dots[0];
    this.previewHolder = preview;
    this.mathMLHolder = mathML;
    this.latexHolder = latex;
    this.ctx = this.canvas.getContext('2d');
    this.strokes = [];
    this.stroke = [];
    this.clicked = 0;
  }

  Drawing.prototype.getMousePos = function(e) {
    var pos, rect;
    rect = this.canvas.getBoundingClientRect();
    return pos = {
      x: e.clientX - rect.left,
      y: e.clientY - rect.top
    };
  };

  Drawing.prototype.start = function(e) {
    var p;
    p = this.getMousePos(e);
    this.clicked = 1;
    this.ctx.beginPath();
    this.ctx.moveTo(p.x, p.y);
    return this.stroke = [[p.x, p.y]];
  };

  Drawing.prototype.move = function(e) {
    var p;
    if (this.clicked) {
      p = this.getMousePos(e);
      this.ctx.lineTo(p.x, p.y);
      this.ctx.stroke();
      return this.stroke.push([p.x, p.y]);
    }
  };

  Drawing.prototype.stop = function(e) {
    this.clicked = 0;
    return this.strokes.push(this.stroke);
  };

  Drawing.prototype.init = function() {
    this.ctx.fillStyle = '#FFFFFF';
    this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    this.ctx.lineWidth = 3;
    this.ctx.lineCap = 'round';
    $(this.canvas).on('mousedown', (function(_this) {
      return function(e) {
        return _this.start(e);
      };
    })(this));
    $(this.canvas).on('mousemove', (function(_this) {
      return function(e) {
        return _this.move(e);
      };
    })(this));
    return $(this.canvas).on('mouseup', (function(_this) {
      return function(e) {
        return _this.stop(e);
      };
    })(this));
  };

  Drawing.prototype.recognize = function() {
    var data, thiz;
    thiz = this;
    $(this.dots).hide();
    $(this.canvas).show();
    data = JSON.stringify(this.strokes);
    return $.ajax({
      type: 'POST',
      contentType: 'application/json; charset=utf-8',
      url: url_recognize,
      data: data,
      dataType: 'json'
    }).done(function(result) {
      var latex, mathML;
      mathML = result[0];
      latex = result[1];
      return thiz.display(mathML, latex);
    }).fail(function(response) {
      return alert(response.responseText);
    });
  };

  Drawing.prototype.display = function(mathML, latex) {
    this.mathMLHolder.text(mathML);
    this.latexHolder.text(latex);
    this.previewHolder.html('$' + latex + '$');
    return MathJax.Hub.Queue(['Typeset', MathJax.Hub, this.previewHolder[0]]);
  };

  Drawing.prototype.drawDots = function(button) {
    var color, ctx, getRandomColor, point, stroke, x, y, _i, _len, _ref, _results;
    $(this.canvas).hide();
    $(this.dots).show();
    ctx = this.dots.getContext('2d');
    ctx.fillStyle = '#EEEEEE';
    ctx.fillRect(0, 0, this.dots.width, this.dots.height);
    getRandomColor = function() {
      var color, i, letters, _i;
      letters = '0123456789ABCDEF'.split('');
      color = '#';
      for (i = _i = 0; _i < 6; i = ++_i) {
        color += letters[Math.floor(Math.random() * 16)];
      }
      return color;
    };
    _ref = this.strokes;
    _results = [];
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      stroke = _ref[_i];
      color = getRandomColor();
      _results.push((function() {
        var _j, _len1, _results1;
        _results1 = [];
        for (_j = 0, _len1 = stroke.length; _j < _len1; _j++) {
          point = stroke[_j];
          x = point[0];
          y = point[1];
          ctx.fillStyle = color;
          _results1.push(ctx.fillRect(x, y, 5, 5));
        }
        return _results1;
      })());
    }
    return _results;
  };

  Drawing.prototype.clean = function() {
    this.clicked = 0;
    this.stroke = [];
    this.strokes = [];
    this.ctx.fillStyle = '#FFFFFF';
    this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    this.previewHolder.html('');
    this.latexHolder.html('');
    return this.mathMLHolder.html('');
  };

  return Drawing;

})();
