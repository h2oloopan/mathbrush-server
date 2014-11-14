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
	

	//initialize
	scg::SetVerbosity(0);
	scg::SetProfilePath("");
	scg::SetUserProfilePath("");
	scg::SetTabletResolution(132);

	if (scg::InitializeRecognizer()) {
		return scope.Close(False());
	}


	scg::SetUserProfile("");
	scg::MathRecognizer* recognizer = scg::CreateMathRecognizer();

	log("* Initialization Done");

	//recognize
	log("* Input:");

	scg::RawStroke* rawStrokes = new scg::RawStroke[strokes->Length()];

	for (int i = 0; i < strokes->Length(); i++) {
		log("* Stroke:");
		Local<Array> stroke = Local<Array>::Cast(strokes->Get(i));

		long* xs = new long[stroke->Length()];
		long* ys = new long[stroke->Length()];

		for (int j = 0; j < stroke->Length(); j++) {
			Local<Array> points = Local<Array>::Cast(stroke->Get(j));
			Local<Number> x = Local<Number>::Cast(points->Get(0));
			Local<Number> y = Local<Number>::Cast(points->Get(1));
			log(std::string("- x: ") + std::string(*(String::Utf8Value(x->ToString()))) + 
				std::string(" y: ") + std::string(*(String::Utf8Value(y->ToString()))));

			xs[j] = (long)x->NumberValue();
			ys[j] = (long)y->NumberValue();
		}

		rawStrokes[i].set_points(xs, ys, stroke->Length());
	}
	recognizer->AddStrokes(rawStrokes, strokes->Length());
	delete[] rawStrokes;

	scg::ExpressionTree* tree = (scg::ExpressionTree*)recognizer->GetTopExpression();
	if (tree == NULL) return scope.Close(False());





	//return mathML
	std::string mathML = "";
	std::string latex = "";
	if (tree->HasLongForm()) {
		scg::ExpressionIterator* iterator = tree->CreateLongFormIterator();
		scg::ExpressionTree* expr = (scg::ExpressionTree*) iterator->next();
		iterator->release();
		if (expr != NULL) {
			mathML = expr->long_str();
			latex = expr->latex_str();
			expr->release();
		}
	}
	else {
		mathML = tree->long_str();
		latex = tree->latex_str();
	}

	log("* MathML: " + mathML);
	log("* Latex: " + latex);

	Local<String> mathMLStr = String::New(mathML.c_str());
	Local<String> latexStr = String::New(latex.c_str());
	Local<Array> result = Array::New(2);
	result->Set(0, mathMLStr);
	result->Set(1, latexStr);
	return scope.Close(result);
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