# cuMES parsing utilities

These files are copied from sibling project `cuMES` at commit `f7036ab` under
its MIT license:

- `clap.h`: command-line argument parser;
- `JsonParser.h`: JSON parser, with its include guard and namespace changed
  to `meow::json` to avoid an ODR collision when meow links cuMES.

Behavioral changes should be made upstream first and synchronized here.

