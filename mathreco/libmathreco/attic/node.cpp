#include "node.h"


namespace scg
{


void
ExpressionNode::AddTag(const std::string &tag)
{
    tags.push_back(tag);
}


bool
ExpressionNode::HasTag(const std::string &tag)
{
    return std::find(tags.begin(), tags.end(), tag) != tags.end();
}


void
ExpressionNode::RemoveTag(const std::string &tag)
{
    std::list<std::string>::iterator i;
    i = std::find(tags.begin(), tags.end(), tag);
    if (i != tags.end()) {
        tags.erase(i);
    }
}



void
ExpressionNode::LinkTo(ExpressionNode *to, const std::string &name, double parm)
{
    links.push_back(DEBUG_NEW ExpressionLink(this, to, name, parm));
}


bool
ExpressionNode::HasLink(ExpressionNode *to, const std::string &name)
{
    for (std::vector<ExpressionLink *>::iterator i = links.begin(); i != links.end(); ++i) {
        if ((!to || (*i)->to == to) && (*i)->name == name) {
            return true;
        }
    }
    return false;
}


void
ExpressionNode::RemoveLink(ExpressionNode *to, const std::string &name)
{
    for (std::vector<ExpressionLink *>::iterator i = links.begin(); i != links.end(); ++i) {
        if ((*i)->to == to && (*i)->name == name) {
            links.erase(i);
            break;
        }
    }
}


LinkCollection
ExpressionNode::GetLinks()
{
    return LinkCollection(links);
}


void
ExpressionNode::AddMatch(const Match &match)
{
    reco_results.push_back(match);
}


}

