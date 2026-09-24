"""
Cross Section Generator Plugin for QGIS
"""

def classFactory(iface):
    from .cross_section_tool import CrossSectionPlugin
    return CrossSectionPlugin(iface)