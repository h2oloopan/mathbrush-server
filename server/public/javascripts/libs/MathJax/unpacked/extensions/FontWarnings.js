!function(e,t){var o="http://www.stixfonts.org/",a="https://github.com/mathjax/MathJax/tree/master/fonts/HTML-CSS/TeX/otf",s=e.CombineConfig("FontWarnings",{messageStyle:{position:"fixed",bottom:"4em",left:"3em",width:"40em",border:"3px solid #880000","background-color":"#E0E0E0",color:"black",padding:"1em","font-size":"small","white-space":"normal","border-radius":".75em","-webkit-border-radius":".75em","-moz-border-radius":".75em","-khtml-border-radius":".75em","box-shadow":"4px 4px 10px #AAAAAA","-webkit-box-shadow":"4px 4px 10px #AAAAAA","-moz-box-shadow":"4px 4px 10px #AAAAAA","-khtml-box-shadow":"4px 4px 10px #AAAAAA",filter:"progid:DXImageTransform.Microsoft.dropshadow(OffX=3, OffY=3, Color='gray', Positive='true')"},Message:{webFont:[["closeBox"],["webFont","MathJax is using web-based fonts to display the mathematics on this page.  These take time to download, so the page would render faster if you installed math fonts directly in your system's font folder."],["fonts"]],imageFonts:[["closeBox"],["imageFonts","MathJax is using its image fonts rather than local or web-based fonts. This will render slower than usual, and the mathematics may not print at the full resolution of your printer."],["fonts"],["webFonts"]],noFonts:[["closeBox"],["noFonts","MathJax is unable to locate a font to use to display its mathematics, and image fonts are not available, so it is falling back on generic unicode characters in hopes that your browser will be able to display them.  Some characters may not show up properly, or possibly not at all."],["fonts"],["webFonts"]]},HTML:{closeBox:[["div",{style:{position:"absolute",overflow:"hidden",top:".1em",right:".1em",border:"1px outset",width:"1em",height:"1em","text-align":"center",cursor:"pointer","background-color":"#EEEEEE",color:"#606060","border-radius":".5em","-webkit-border-radius":".5em","-moz-border-radius":".5em","-khtml-border-radius":".5em"},onclick:function(){n.div&&0===n.fade&&(n.timer&&clearTimeout(n.timer),n.div.style.display="none")}},[["span",{style:{position:"relative",bottom:".2em"}},["x"]]]]],webFonts:[["p"],["webFonts","Most modern browsers allow for fonts to be downloaded over the web. Updating to a more recent version of your browser (or changing browsers) could improve the quality of the mathematics on this page."]],fonts:[["p"],["fonts","MathJax can use either the [STIX fonts](%1) or the [MathJax TeX fonts](%2).  Download and install one of those fonts to improve your MathJax experience.",o,a]],STIXfonts:[["p"],["STIXPage","This page is designed to use the [STIX fonts](%1).  Download and install those fonts to improve your MathJax experience.",o]],TeXfonts:[["p"],["TeXPage","This page is designed to use the [MathJax TeX fonts](%1).  Download and install those fonts to improve your MathJax experience.",a]]},removeAfter:12e3,fadeoutSteps:10,fadeoutTime:1500});MathJax.Hub.Browser.isIE9&&document.documentMode>=9&&delete s.messageStyle.filter;var n={div:null,fade:0},i=function(o){if(!n.div){var a=MathJax.OutputJax["HTML-CSS"],i=document.body;e.Browser.isMSIE?"fixed"===s.messageStyle.position&&(MathJax.Message.Init(),i=document.getElementById("MathJax_MSIE_Frame")||i,i!==document.body&&(s.messageStyle.position="absolute")):delete s.messageStyle.filter,s.messageStyle.maxWidth=document.body.clientWidth-75+"px";for(var l=0;l<o.length;)if(o[l]instanceof Array)if(1===o[l].length&&s.HTML[o[l][0]])o.splice.apply(o,[l,1].concat(s.HTML[o[l][0]]));else if("string"==typeof o[l][1]){var d=MathJax.Localization.lookupPhrase(["FontWarnings",o[l][0]],o[l][1]);d=MathJax.Localization.processMarkdown(d,o[l].slice(2),"FontWarnings"),o.splice.apply(o,[l,1].concat(d)),l+=d.length}else l++;else l++;n.div=a.addElement(i,"div",{id:"MathJax_FontWarning",style:s.messageStyle},o),MathJax.Localization.setCSS(n.div),s.removeAfter&&e.Register.StartupHook("End",function(){n.timer=setTimeout(r,s.removeAfter)}),t.Cookie.Set("fontWarn",{warned:!0})}},r=function(){if(n.fade++,n.timer&&delete n.timer,n.fade<s.fadeoutSteps){var e=1-n.fade/s.fadeoutSteps;n.div.style.opacity=e,n.div.style.filter="alpha(opacity="+Math.floor(100*e)+")",setTimeout(r,s.fadeoutTime/s.fadeoutSteps)}else n.div.style.display="none"};t.Cookie.Get("fontWarn").warned||e.Startup.signal.Interest(function(e){if(e.match(/HTML-CSS Jax - /)&&!n.div){var t,o=MathJax.OutputJax["HTML-CSS"],a=o.config.availableFonts,r=a&&a.length;r?1===a.length&&(s.HTML.fonts=s.HTML[a[0]+"fonts"]):s.HTML.fonts=[""],o.allowWebFonts&&(s.HTML.webfonts=[""]),e.match(/- Web-Font/)?r&&(t="webFont"):e.match(/- using image fonts/)?t="imageFonts":e.match(/- no valid font/)&&(t="noFonts"),t&&s.Message[t]&&MathJax.Localization.loadDomain("FontWarnings",[i,s.Message[t]])}})}(MathJax.Hub,MathJax.HTML),MathJax.Ajax.loadComplete("[MathJax]/extensions/FontWarnings.js");