"""Exceptions for asap_to_kite."""


class AsapToKiteError(Exception):
    """Base exception."""


class InvalidInputError(AsapToKiteError):
    """Invalid or missing input files."""


class ProcessingError(AsapToKiteError):
    """Processing failed."""


class ValidationError(AsapToKiteError):
    """Validation failed."""
