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
			"libraries": ["../../../mathreco/lib/libmathreco.a"]
		}
	]
}