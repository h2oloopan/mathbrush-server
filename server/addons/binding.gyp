{
	"targets": [
		{
			"target_name": "hello",
			"sources": [ "hello.cc" ]
		},
		{
			"target_name": "recognizer",
			"sources": [ "recognizer.cc" ],
			"include_dirs": ["../../mathreco/include"],
			"conditions": [
				"OS=='mac'",
				{
					"libraries": ["../../../libs/libmathreco_osx.a"]
				}
			]
		}
	]
}