from shiny import App
from fengalotl.mod.ui import app_ui as _ui
from fengalotl.mod.server import server as _server

app = App(_ui, _server)
