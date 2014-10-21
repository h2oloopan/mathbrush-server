#define BUILDING_NODE_EXTENSION

#include <node.h>
#include <v8.h>
#include <cstring>
#include <iostream>

#include "MathRecognizer.h"

using namespace v8;

Persistent<Object> logger;

void log(const std::string msg) {
	const unsigned argc = 1;
	Local<Value> argv[argc] = {
		Local<Value>::New(String::New(msg.c_str()))
	};
	Local<Function> func = Local<Function>::Cast(logger->Get(String::NewSymbol("log")));
	func->Call(Context::GetCurrent()->Global(), argc, argv);
}

void log(const char* msg) {
	const unsigned argc = 1;
	Local<Value> argv[argc] = {
		Local<Value>::New(String::New(msg))
	};
	Local<Function> func = Local<Function>::Cast(logger->Get(String::NewSymbol("log")));
	func->Call(Context::GetCurrent()->Global(), argc, argv);
}

//input signature: name, strokes, return latex
Handle<Value> Recognize(const Arguments& args) {
	HandleScope scope;
	Local<String> name = String::New(*String::Utf8Value(args[0]->ToString()));
	Local<Array> strokes = Local<Array>::Cast(args[1]);
	String::Utf8Value str(name);
	log("*********************");
	log("* Start Recognition");
	log(std::string("* User: ") + std::string(*str));
	log("* Input:");
	for (int i = 0; i < strokes->Length(); i++) {
		log("* Stroke:");
		Local<Array> stroke = Local<Array>::Cast(strokes->Get(i));
		for (int j = 0; j < stroke->Length(); j++) {
			Local<Array> points = Local<Array>::Cast(stroke->Get(j));
			Local<Number> x = Local<Number>::Cast(points->Get(0));
			Local<Number> y = Local<Number>::Cast(points->Get(1));
			log(std::string("- x: ") + std::string(*(String::Utf8Value(x->ToString()))) + 
				std::string(" y: ") + std::string(*(String::Utf8Value(y->ToString()))));
		}
	}

	return scope.Close(name);
}

Handle<Value> Initialize(const Arguments& args) {
	HandleScope scope;
	logger = Persistent<Object>::New(args[0]->ToObject());
	return scope.Close(Null());
}

void Init(Handle<Object> exports) {
	exports->Set(String::NewSymbol("recognize"), FunctionTemplate::New(Recognize)->GetFunction());
	exports->Set(String::NewSymbol("initialize"), FunctionTemplate::New(Initialize)->GetFunction());
}

NODE_MODULE(recognizer, Init)