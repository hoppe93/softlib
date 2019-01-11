/**
 * Implementation of the configuration interpreter
 */

#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>
#include <softlib/Configuration.h>
#include <softlib/ConfigurationException.h>

using namespace std;

/**
 * Moves one step ahead in the token stream.
 * Returns TRUE as long as the end of the stream
 * (i.e. EOF) has not been reached.
 */
bool Configuration::advance() {
	// If current token is EOF, we're already done
	if (tokent() == CTKN_EOF) return false;

	ConfigToken *tkn = *(++this->tknstream);
	return (tkn->GetType()!=CTKN_EOF);
}
/**
 * Advances the stream until the desired token type
 * lies ahead, or until the next token is EOF. Stops
 * right before the next token, so that the next token
 * is guaranteed to be of one of the types passed to it if
 * this function returns TRUE, and EOF if it returns
 * FALSE. Any CTKN_INVALID is ignored.
 *
 * type: Array of types to search for. The function stops
 *   at the first match in this array.
 */
bool Configuration::advanceTo(std::initializer_list<conftoken_t> args) {
	do {
		ConfigToken *tkn = *(this->tknstream+1);

		for (auto type : args) {
			if (tkn->GetType() == type)
				return true;
		}
	} while (advance());

	return false;
}

/**
 * Returns the current character position.
 */
size_t Configuration::currcharpos() {
    return (token()->GetCharpos());
}
/**
 * Returns the name of the current file.
 */
string Configuration::currfile() {
    return (token()->GetFilename());
}
/**
 * Returns the current line number.
 */
size_t Configuration::currline() {
    return (token()->GetLine());
}

/**
 * Advance to next token in the stream, and quit
 * with an error message if it is not of type 'type'.
 *
 * type: Expected type of the next token.
 */
void Configuration::expect(conftoken_t type) {
	if (!advance() || tokent() != type) {
		const string tps  = ttos(type);
		const string toks = ttos(tokent());
		throw ConfigurationException(currline(), currcharpos(), currfile(), "Expected " + tps + " but got " + toks);
	}
}
/**
 * Return the current token in the stream
 */
ConfigToken *Configuration::token() {
	ConfigToken *tkn = *(this->tknstream);
	return tkn;
}
/**
 * Return the type of the current token
 */
conftoken_t Configuration::tokent() {
	ConfigToken *tkn = *(this->tknstream);
	return tkn->GetType();
}

/**
 * Return the type of the next token
 */
conftoken_t Configuration::tokent_n() {
	// If EOF has been reached, next token is EOF
	if (tokent() == CTKN_EOF) return CTKN_EOF;

	ConfigToken *tkn = *(this->tknstream+1);
	return tkn->GetType();
}

/** Convert a token-type to a string.
 *
 * type: The token-type to translate to a string.
 */
const char *Configuration::ttos(conftoken_t type) {
	switch (type) {
		case CTKN_NAME: return "name";
		case CTKN_VALUE: return "value";
		case CTKN_EQUALS: return "=";
        case CTKN_LPAR: return "(";
        case CTKN_RPAR: return ")";
		case CTKN_BLOCKTYPE: return "block-type";
		case CTKN_BLOCKNAME: return "block-name";
		case CTKN_BLOCK_START: return "{";
		case CTKN_BLOCK_END: return "}";
		case CTKN_END_STATEMENT: return ";";
		case CTKN_EOF: return "<EOF>";
		default: return "<invalid>";
	}
}

/**
 * Convert the stream of tokens to a 'ConfigBlock'.
 */
ConfigBlock Configuration::Interpret(vector<ConfigToken*>& tkns) {
	ConfigBlock rootblock = ConfigBlock(0, "<root>");
	ConfigBlock *blck;
	ConfigToken *tkn;
	confblock_t blockt;
	string name, value, sndtype;
	this->tknstream = tkns.begin();

	blck = &rootblock;

	do {
		try {
			switch (tokent()) {
				case CTKN_NAME: {
					tkn = token();
					name = tkn->GetValue();

					expect(CTKN_EQUALS);
					expect(CTKN_VALUE);
					tkn = token();
					value = tkn->GetValue();
					Setting& set = blck->AddSetting(name, value);

					// Append any trailing values to build a vector
					while (tokent_n() == CTKN_VALUE) {
						advance();
						tkn = token();
						value = tkn->GetValue();
						set.AppendValue(value);
					}

					expect(CTKN_END_STATEMENT);
				} break;
				case CTKN_BLOCKTYPE: {
					tkn = token();
					blockt = StringToBlockType(tkn->GetValue());

					expect(CTKN_BLOCKNAME);
					tkn = token();
					name = tkn->GetValue();

                    if (tokent_n() == CTKN_LPAR) {
                        advance();
                        expect(CTKN_NAME);
                        tkn = token();
                        sndtype = tkn->GetValue();
                        expect(CTKN_RPAR);
                    }

					expect(CTKN_BLOCK_START);
					
					// Add the new sub-block and make it the active block
					blck = blck->AddSubBlock(blockt, name, sndtype);
				} break;
				case CTKN_BLOCK_END: {
					if (blck->HasParent())
						blck = blck->GetParent();
					else throw ConfigurationException(currline(), currcharpos(), currfile(), "Syntax error: Mismatched '{'.");
				} break;
				default:
					const string toks = ttos(tokent());
					throw ConfigurationException(currline(), currcharpos(), currfile(), "Syntax error: Expected name or block-type, but got " + toks);
			}
		} catch (ConfigurationException& ex) {
			// Advance to next statement
			cerr << "ERROR: " << ex.whats() << endl;
			advanceTo({CTKN_NAME, CTKN_BLOCKTYPE, CTKN_BLOCK_END});
            this->SetErrorFlag();
		}
	} while (advance());

	// Have we returned to the root block?
	if (blck != &rootblock) {
		throw ConfigurationException(currline(), currcharpos(), currfile(), "Mismatched '{'. Document was not closed.");
	}

	return rootblock;
}

