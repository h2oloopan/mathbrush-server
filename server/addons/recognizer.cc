#define BUILDING_NODE_EXTENSION

#include <node.h>
#include <v8.h>
#include <cstring>
#include <iostream>

#include "MathRecognizer.h"

using namespace v8;

static Local<Function> logger;

void log(const std::string msg) {
	const unsigned argc = 1;
	Local<Value> argv[argc] = {
		Local<Value>::New(String::New(msg.c_str()))
	};
	logger->Call(Context::GetCurrent()->Global(), argc, argv);
}

//input signature: name, strokes, return latex
Handle<Value> Recognize(const Arguments& args) {
	HandleScope scope;
	Local<String> name = String::New(*String::Utf8Value(args[0]->ToString()));
	Local<Array> strokes = Local<Array>::Cast(args[1]);
	String::Utf8Value str(name);
	std::cout << std::string(*str);
	return scope.Close(name);
}

Handle<Value> Initialize(const Arguments& args) {
	HandleScope scope;
	logger = Local<Function>::Cast(args[0]);
}

void Init(Handle<Object> exports) {
	exports->Set(String::NewSymbol("recognize"), FunctionTemplate::New(Recognize)->GetFunction());
	exports->Set(String::NewSymbol("initialize"), FunctionTemplate::New(Initialize)->GetFunction());
}

NODE_MODULE(recognizer, Init)