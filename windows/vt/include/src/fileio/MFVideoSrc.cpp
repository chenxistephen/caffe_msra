//+----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class to read in video frames using Media Foundation (MF).
//
//  History:
//      2010/06/25-t-parasa
//          Created
//      2011/07/17-sbaker
//          Added support for Video Images C*VImg
//          Added CMFVideoDst
//
//-----------------------------------------------------------------------------

#include "stdafx.h"

using namespace vt;

// Misc
#include <assert.h>
#include <memory>

#include "vt_io.h"
#include "MFVideoSrc.h"

#if ( defined(ENABLE_HARDWARE_DECODER) && (_MSC_VER >= 1700) )
#pragma comment (lib, "d3d11")
#endif

#if defined(VT_WINRT)
using namespace Windows::Storage;
using namespace Windows::Storage::Streams;
#endif

#pragma warning ( disable : 6326 ) // comparing two constants
// MF error codes
struct MF_ERROR_MAP
{
    LPCWSTR wszMFError;
    HRESULT hr;
};

bool vt::GetMFErrorString(HRESULT hr, wstring& errmsg)
{
    const static MF_ERROR_MAP MFErrorMap[] =
    {
        {L"Platform not initialized. Please call MFStartup().", MF_E_PLATFORM_NOT_INITIALIZED},
        {L"The buffer was too small to carry out the requested action.", MF_E_BUFFERTOOSMALL},
        {L"The request is invalid in the current state.", MF_E_INVALIDREQUEST},
        {L"The stream number provided was invalid.", MF_E_INVALIDSTREAMNUMBER},
        {L"The data specified for the media type is invalid, inconsistent, or not supported by this object.", MF_E_INVALIDMEDIATYPE},
        {L"The callee is currently not accepting further input.", MF_E_NOTACCEPTING},
        {L"This object needs to be initialized before the requested operation can be carried out.", MF_E_NOT_INITIALIZED},
        {L"The requested representation is not supported by this object.", MF_E_UNSUPPORTED_REPRESENTATION},
        {L"An object ran out of media types to suggest therefore the requested chain of streaming objects cannot be completed.", MF_E_NO_MORE_TYPES},
        {L"The object does not support the specified service.", MF_E_UNSUPPORTED_SERVICE},
        {L"An unexpected error has occurred in the operation requested.", MF_E_UNEXPECTED},
        {L"Invalid name.", MF_E_INVALIDNAME},
        {L"Invalid type.", MF_E_INVALIDTYPE},
        {L"Invalid type.", MF_E_INVALID_FILE_FORMAT},
        {L"Invalid index.", MF_E_INVALIDINDEX},
        {L"An invalid timestamp was given.", MF_E_INVALID_TIMESTAMP},
        {L"The scheme of the given URL is unsupported.", MF_E_UNSUPPORTED_SCHEME},
        {L"The byte stream type of the given URL is unsupported.", MF_E_UNSUPPORTED_BYTESTREAM_TYPE},
        {L"The given time format is unsupported.", MF_E_UNSUPPORTED_TIME_FORMAT},
        {L"The Media Sample does not have a timestamp.", MF_E_NO_SAMPLE_TIMESTAMP},
        {L"The Media Sample does not have a duration.", MF_E_NO_SAMPLE_DURATION},
        {L"The request failed because the data in the stream is corrupt.\n.", MF_E_INVALID_STREAM_DATA},
        {L"Real time services are not available.", MF_E_RT_UNAVAILABLE},
        {L"The specified rate is not supported.", MF_E_UNSUPPORTED_RATE},
        {L"This component does not support stream-thinning.", MF_E_THINNING_UNSUPPORTED},
        {L"The call failed because no reverse playback rates are available.", MF_E_REVERSE_UNSUPPORTED},
        {L"The requested rate transition cannot occur in the current state.", MF_E_UNSUPPORTED_RATE_TRANSITION},
        {L"The requested rate change has been pre-empted and will not occur.", MF_E_RATE_CHANGE_PREEMPTED},
        {L"The specified object or value does not exist.", MF_E_NOT_FOUND},
        {L"The requested value is not available.", MF_E_NOT_AVAILABLE},
        {L"The specified operation requires a clock and no clock is available.", MF_E_NO_CLOCK},
        {L"This callback and state had already been passed in to this event generator earlier.", MF_S_MULTIPLE_BEGIN},
        {L"This callback has already been passed in to this event generator.", MF_E_MULTIPLE_BEGIN},
        {L"Some component is already listening to events on this event generator.", MF_E_MULTIPLE_SUBSCRIBERS},
        {L"This timer was orphaned before its callback time arrived.", MF_E_TIMER_ORPHANED},
        {L"A state transition is already pending.", MF_E_STATE_TRANSITION_PENDING},
        {L"The requested state transition is unsupported.", MF_E_UNSUPPORTED_STATE_TRANSITION},
        {L"An unrecoverable error has occurred.", MF_E_UNRECOVERABLE_ERROR_OCCURRED},
        {L"The provided sample has too many buffers.", MF_E_SAMPLE_HAS_TOO_MANY_BUFFERS},
        {L"The provided sample is not writable.", MF_E_SAMPLE_NOT_WRITABLE},
        {L"The provided sample is not writable.", MF_E_INVALID_KEY},
        {L"The provided sample is not writable.", MF_E_BAD_STARTUP_VERSION},
        {L"The caption of the given URL is unsupported.", MF_E_UNSUPPORTED_CAPTION},
        {L"The operation on the current offset is not permitted.", MF_E_INVALID_POSITION},
        {L"The requested attribute was not found.", MF_E_ATTRIBUTENOTFOUND},
        {L"The specified property type is not allowed in this context.", MF_E_PROPERTY_TYPE_NOT_ALLOWED},
        {L"The specified property type is not supported.", MF_E_PROPERTY_TYPE_NOT_SUPPORTED},
        {L"The specified property is empty.", MF_E_PROPERTY_EMPTY},
        {L"The specified property is not empty.", MF_E_PROPERTY_NOT_EMPTY},
        {L"The vector property specified is not allowed in this context.", MF_E_PROPERTY_VECTOR_NOT_ALLOWED},
        {L"A vector property is required in this context.", MF_E_PROPERTY_VECTOR_REQUIRED},
        {L"The operation is cancelled.", MF_E_OPERATION_CANCELLED},
        {L"The provided bytestream was expected to be seekable and it is not.", MF_E_BYTESTREAM_NOT_SEEKABLE},
        {L"The Media Foundation platform is disabled when the system is running in Safe Mode.", MF_E_DISABLED_IN_SAFEMODE},
        {L"The Media Source could not parse the byte stream.", MF_E_CANNOT_PARSE_BYTESTREAM},
        {L"Mutually exclusive flags have been specified to source resolver. This flag combination is invalid.", MF_E_SOURCERESOLVER_MUTUALLY_EXCLUSIVE_FLAGS},
        {L"MediaProc is in the wrong state", MF_E_MEDIAPROC_WRONGSTATE},
        {L"Real time I/O service can not provide requested throughput.", MF_E_RT_THROUGHPUT_NOT_AVAILABLE},
        {L"The workqueue cannot be registered with more classes.", MF_E_RT_TOO_MANY_CLASSES},
        {L"This operation cannot succeed because another thread owns this object.", MF_E_RT_WOULDBLOCK},
        {L"Internal. Bitpump not found.", MF_E_NO_BITPUMP},
        {L"No more RT memory available.", MF_E_RT_OUTOFMEMORY},
        {L"An MMCSS class has not been set for this work queue.", MF_E_RT_WORKQUEUE_CLASS_NOT_SPECIFIED},
        {L"Insufficient memory for response.", MF_E_INSUFFICIENT_BUFFER},
        {L"Activate failed to create mediasink. Call OutputNode::GetUINT32(MF_TOPONODE_MAJORTYPE) for more information. ", MF_E_CANNOT_CREATE_SINK},
        {L"The length of the provided bytestream is unknown.", MF_E_BYTESTREAM_UNKNOWN_LENGTH},
        {L"The media session cannot pause from a stopped state.", MF_E_SESSION_PAUSEWHILESTOPPED},
        {L"The activate could not be created in the remote process for some reason it was replaced with empty one.", MF_S_ACTIVATE_REPLACED},
        {L"The data specified for the media type is supported, but would require a format change, which is not supported by this object.", MF_E_FORMAT_CHANGE_NOT_SUPPORTED},
        {L"The operation failed because an invalid combination of workqueue ID and flags was specified.", MF_E_INVALID_WORKQUEUE},
        {L"No DRM support is available.", MF_E_DRM_UNSUPPORTED},
        {L"This operation is not authorized.", MF_E_UNAUTHORIZED},
        {L"The value is not in the specified or valid range.", MF_E_OUT_OF_RANGE},
        {L"The registered codec merit is not valid.", MF_E_INVALID_CODEC_MERIT},
        {L"Hardware MFT failed to start streaming due to lack of hardware resources.", MF_E_HW_MFT_FAILED_START_STREAMING},
        {L"Hardware MFT failed to start streaming due to lack of hardware resources.", MF_E_HW_MFT_FAILED_START_STREAMING},
        {L"Hardware MFT failed to start streaming due to lack of hardware resources.", MF_E_HW_MFT_FAILED_START_STREAMING},
        {L"Parsing is still in progress and is not yet complete.", MF_S_ASF_PARSEINPROGRESS},
        {L"Parsing is still in progress and is not yet complete.", MF_S_ASF_PARSEINPROGRESS},
        {L"Parsing is still in progress and is not yet complete.", MF_S_ASF_PARSEINPROGRESS},
        {L"Not enough data have been parsed to carry out the requested action.", MF_E_ASF_PARSINGINCOMPLETE},
        {L"There is a gap in the ASF data provided.", MF_E_ASF_MISSINGDATA},
        {L"The data provided are not valid ASF.", MF_E_ASF_INVALIDDATA},
        {L"The packet is opaque, so the requested information cannot be returned.", MF_E_ASF_OPAQUEPACKET},
        {L"The requested operation failed since there is no appropriate ASF index.", MF_E_ASF_NOINDEX},
        {L"The value supplied is out of range for this operation.", MF_E_ASF_OUTOFRANGE},
        {L"The index entry requested needs to be loaded before it can be available.", MF_E_ASF_INDEXNOTLOADED},
        {L"The packet has reached the maximum number of payloads.", MF_E_ASF_TOO_MANY_PAYLOADS},
        {L"Stream type is not supported.", MF_E_ASF_UNSUPPORTED_STREAM_TYPE},
        {L"One or more ASF packets were dropped.", MF_E_ASF_DROPPED_PACKET},
        {L"One or more ASF packets were dropped.", MF_E_ASF_DROPPED_PACKET},
        {L"One or more ASF packets were dropped.", MF_E_ASF_DROPPED_PACKET},
        {L"There are no events available in the queue.", MF_E_NO_EVENTS_AVAILABLE},
        {L"A media source cannot go from the stopped state to the paused state.", MF_E_INVALID_STATE_TRANSITION},
        {L"The media stream cannot process any more samples because there are no more samples in the stream.", MF_E_END_OF_STREAM},
        {L"The request is invalid because Shutdown() has been called.", MF_E_SHUTDOWN},
        {L"The MP3 object was not found.", MF_E_MP3_NOTFOUND},
        {L"The MP3 parser ran out of data before finding the MP3 object.", MF_E_MP3_OUTOFDATA},
        {L"The file is not really a MP3 file.", MF_E_MP3_NOTMP3},
        {L"The MP3 file is not supported.", MF_E_MP3_NOTSUPPORTED},
        {L"The Media stream has no duration.", MF_E_NO_DURATION},
        {L"The Media format is recognized but is invalid.", MF_E_INVALID_FORMAT},
        {L"The property requested was not found.", MF_E_PROPERTY_NOT_FOUND},
        {L"The property is read only.", MF_E_PROPERTY_READ_ONLY},
        {L"The specified property is not allowed in this context.", MF_E_PROPERTY_NOT_ALLOWED},
        {L"The media source is not started.", MF_E_MEDIA_SOURCE_NOT_STARTED},
        {L"The Media format is recognized but not supported.", MF_E_UNSUPPORTED_FORMAT},
        {L"The MPEG frame has bad CRC.", MF_E_MP3_BAD_CRC},
        {L"The file is not protected.", MF_E_NOT_PROTECTED},
        {L"The media source is in the wrong state", MF_E_MEDIA_SOURCE_WRONGSTATE},
        {L"No streams are selected in source presentation descriptor.", MF_E_MEDIA_SOURCE_NO_STREAMS_SELECTED},
        {L"No key frame sample was found.", MF_E_CANNOT_FIND_KEYFRAME_SAMPLE},
        {L"No key frame sample was found.", MF_E_CANNOT_FIND_KEYFRAME_SAMPLE},
        {L"No key frame sample was found.", MF_E_CANNOT_FIND_KEYFRAME_SAMPLE},
        {L"An attempt to acquire a network resource failed.", MF_E_NETWORK_RESOURCE_FAILURE},
        {L"Error writing to the network.", MF_E_NET_WRITE},
        {L"Error reading from the network.", MF_E_NET_READ},
        {L"Internal. Entry cannot complete operation without network.", MF_E_NET_REQUIRE_NETWORK},
        {L"Internal. Async op is required.", MF_E_NET_REQUIRE_ASYNC},
        {L"Internal. Bandwidth levels are not supported.", MF_E_NET_BWLEVEL_NOT_SUPPORTED},
        {L"Internal. Stream groups are not supported.", MF_E_NET_STREAMGROUPS_NOT_SUPPORTED},
        {L"Manual stream selection is not supported.", MF_E_NET_MANUALSS_NOT_SUPPORTED},
        {L"Invalid presentation descriptor.", MF_E_NET_INVALID_PRESENTATION_DESCRIPTOR},
        {L"Cannot find cache stream.", MF_E_NET_CACHESTREAM_NOT_FOUND},
        {L"The proxy setting is manual.", MF_I_MANUAL_PROXY},
        {L"Internal. Entry cannot complete operation without input.", MF_E_NET_REQUIRE_INPUT},
        {L"The client redirected to another server.", MF_E_NET_REDIRECT},
        {L"The client is redirected to a proxy server.", MF_E_NET_REDIRECT_TO_PROXY},
        {L"The client reached maximum redirection limit.", MF_E_NET_TOO_MANY_REDIRECTS},
        {L"The server, a computer set up to offer multimedia content to other computers, could not handle your request for multimedia content in a timely manner.  Please try again later.", MF_E_NET_TIMEOUT},
        {L"The control socket is closed by the client.", MF_E_NET_CLIENT_CLOSE},
        {L"The server received invalid data from the client on the control connection.", MF_E_NET_BAD_CONTROL_DATA},
        {L"The server is not a compatible streaming media server.", MF_E_NET_INCOMPATIBLE_SERVER},
        {L"Url.", MF_E_NET_UNSAFE_URL},
        {L"Data is not available.", MF_E_NET_CACHE_NO_DATA},
        {L"End of line.", MF_E_NET_EOL},
        {L"The request could not be understood by the server.", MF_E_NET_BAD_REQUEST},
        {L"The server encountered an unexpected condition which prevented it from fulfilling the request.", MF_E_NET_INTERNAL_SERVER_ERROR},
        {L"Session not found.", MF_E_NET_SESSION_NOT_FOUND},
        {L"There is no connection established with the Windows Media server. The operation failed.", MF_E_NET_NOCONNECTION},
        {L"The network connection has failed.", MF_E_NET_CONNECTION_FAILURE},
        {L"The Server service that received the HTTP push request is not a compatible version of Windows Media Services (WMS).  This error may indicate the push request was received by IIS instead of WMS.  Ensure WMS is started and has the HTTP Server control protocol properly enabled and try again.", MF_E_NET_INCOMPATIBLE_PUSHSERVER},
        {L"The Windows Media server is denying access.  The username and/or password might be incorrect.", MF_E_NET_SERVER_ACCESSDENIED},
        {L"The proxy server is denying access.  The username and/or password might be incorrect.", MF_E_NET_PROXY_ACCESSDENIED},
        {L"Unable to establish a connection to the server.", MF_E_NET_CANNOTCONNECT},
        {L"The specified push template is invalid.", MF_E_NET_INVALID_PUSH_TEMPLATE},
        {L"The specified push publishing point is invalid.", MF_E_NET_INVALID_PUSH_PUBLISHING_POINT},
        {L"The requested resource is in use.", MF_E_NET_BUSY},
        {L"The Publishing Point or file on the Windows Media Server is no longer available.", MF_E_NET_RESOURCE_GONE},
        {L"The proxy experienced an error while attempting to contact the media server.", MF_E_NET_ERROR_FROM_PROXY},
        {L"The proxy did not receive a timely response while attempting to contact the media server.", MF_E_NET_PROXY_TIMEOUT},
        {L"The server is currently unable to handle the request due to a temporary overloading or maintenance of the server.", MF_E_NET_SERVER_UNAVAILABLE},
        {L"The encoding process was unable to keep up with the amount of supplied data.", MF_E_NET_TOO_MUCH_DATA},
        {L"Session not found.", MF_E_NET_SESSION_INVALID},
        {L"The requested URL is not available in offline mode.", MF_E_OFFLINE_MODE},
        {L"A device in the network is blocking UDP traffic.", MF_E_NET_UDP_BLOCKED},
        {L"The specified configuration value is not supported.", MF_E_NET_UNSUPPORTED_CONFIGURATION},
        {L"The networking protocol is disabled.", MF_E_NET_PROTOCOL_DISABLED},
        {L"The networking protocol is disabled.", MF_E_NET_PROTOCOL_DISABLED},
        {L"The networking protocol is disabled.", MF_E_NET_PROTOCOL_DISABLED},
        {L"This object has already been initialized and cannot be re-initialized at this time.", MF_E_ALREADY_INITIALIZED},
        {L"The amount of data passed in exceeds the given bitrate and buffer window.", MF_E_BANDWIDTH_OVERRUN},
        {L"The sample was passed in too late to be correctly processed.", MF_E_LATE_SAMPLE},
        {L"The requested action cannot be carried out until the object is flushed and the queue is emptied.", MF_E_FLUSH_NEEDED},
        {L"The profile is invalid.", MF_E_INVALID_PROFILE},
        {L"The index that is being generated needs to be committed before the requested action can be carried out.", MF_E_INDEX_NOT_COMMITTED},
        {L"The index that is necessary for the requested action is not found.", MF_E_NO_INDEX},
        {L"The requested index cannot be added in-place to the specified ASF content.", MF_E_CANNOT_INDEX_IN_PLACE},
        {L"The ASF leaky bucket parameters must be specified in order to carry out this request.", MF_E_MISSING_ASF_LEAKYBUCKET},
        {L"The stream id is invalid. The valid range for ASF stream id is from 1 to 127.", MF_E_INVALID_ASF_STREAMID},
        {L"The stream id is invalid. The valid range for ASF stream id is from 1 to 127.", MF_E_INVALID_ASF_STREAMID},
        {L"The stream id is invalid. The valid range for ASF stream id is from 1 to 127.", MF_E_INVALID_ASF_STREAMID},
        {L"The requested Stream Sink has been removed and cannot be used.", MF_E_STREAMSINK_REMOVED},
        {L"The various Stream Sinks in this Media Sink are too far out of sync for the requested action to take place.", MF_E_STREAMSINKS_OUT_OF_SYNC},
        {L"Stream Sinks cannot be added to or removed from this Media Sink because its set of streams is fixed.", MF_E_STREAMSINKS_FIXED},
        {L"The given Stream Sink already exists.", MF_E_STREAMSINK_EXISTS},
        {L"Sample allocations have been canceled.", MF_E_SAMPLEALLOCATOR_CANCELED},
        {L"The sample allocator is currently empty, due to outstanding requests.", MF_E_SAMPLEALLOCATOR_EMPTY},
        {L"When we try to sopt a stream sink, it is already stopped ", MF_E_SINK_ALREADYSTOPPED},
        {L"The ASF file sink could not reserve AVIO because the bitrate is unknown.", MF_E_ASF_FILESINK_BITRATE_UNKNOWN},
        {L"No streams are selected in sink presentation descriptor.", MF_E_SINK_NO_STREAMS},
        {L"The sink has not been finalized before shut down. This may cause sink generate a corrupted content.", MF_S_SINK_NOT_FINALIZED},
        {L"A metadata item was too long to write to the output container.", MF_E_METADATA_TOO_LONG},
        {L"The operation failed because no samples were processed by the sink.", MF_E_SINK_NO_SAMPLES_PROCESSED},
        {L"The operation failed because no samples were processed by the sink.", MF_E_SINK_NO_SAMPLES_PROCESSED},
        {L"The operation failed because no samples were processed by the sink.", MF_E_SINK_NO_SAMPLES_PROCESSED},
        {L"There is no available procamp hardware with which to perform color correction.", MF_E_VIDEO_REN_NO_PROCAMP_HW},
        {L"There is no available deinterlacing hardware with which to deinterlace the video stream.", MF_E_VIDEO_REN_NO_DEINTERLACE_HW},
        {L"A video stream requires copy protection to be enabled, but there was a failure in attempting to enable copy protection.", MF_E_VIDEO_REN_COPYPROT_FAILED},
        {L"A component is attempting to access a surface for sharing that is not shared.", MF_E_VIDEO_REN_SURFACE_NOT_SHARED},
        {L"A component is attempting to access a shared device that is already locked by another component.", MF_E_VIDEO_DEVICE_LOCKED},
        {L"The device is no longer available. The handle should be closed and a new one opened.", MF_E_NEW_VIDEO_DEVICE},
        {L"A video sample is not currently queued on a stream that is required for mixing.", MF_E_NO_VIDEO_SAMPLE_AVAILABLE},
        {L"No audio playback device was found.", MF_E_NO_AUDIO_PLAYBACK_DEVICE},
        {L"The requested audio playback device is currently in use.", MF_E_AUDIO_PLAYBACK_DEVICE_IN_USE},
        {L"The audio playback device is no longer present.", MF_E_AUDIO_PLAYBACK_DEVICE_INVALIDATED},
        {L"The audio service is not running.", MF_E_AUDIO_SERVICE_NOT_RUNNING},
        {L"The audio service is not running.", MF_E_AUDIO_SERVICE_NOT_RUNNING},
        {L"The audio service is not running.", MF_E_AUDIO_SERVICE_NOT_RUNNING},
        {L"The topology contains an invalid optional node.  Possible reasons are incorrect number of outputs and inputs or optional node is at the beginning or end of a segment. ", MF_E_TOPO_INVALID_OPTIONAL_NODE},
        {L"No suitable transform was found to decrypt the content. ", MF_E_TOPO_CANNOT_FIND_DECRYPTOR},
        {L"No suitable transform was found to encode or decode the content. ", MF_E_TOPO_CODEC_NOT_FOUND},
        {L"Unable to find a way to connect nodes", MF_E_TOPO_CANNOT_CONNECT},
        {L"Unsupported operations in topoloader", MF_E_TOPO_UNSUPPORTED},
        {L"The topology or its nodes contain incorrectly set time attributes", MF_E_TOPO_INVALID_TIME_ATTRIBUTES},
        {L"The topology contains loops, which are unsupported in media foundation topologies", MF_E_TOPO_LOOPS_IN_TOPOLOGY},
        {L"A source stream node in the topology does not have a presentation descriptor", MF_E_TOPO_MISSING_PRESENTATION_DESCRIPTOR},
        {L"A source stream node in the topology does not have a stream descriptor", MF_E_TOPO_MISSING_STREAM_DESCRIPTOR},
        {L"A stream descriptor was set on a source stream node but it was not selected on the presentation descriptor", MF_E_TOPO_STREAM_DESCRIPTOR_NOT_SELECTED},
        {L"A source stream node in the topology does not have a source", MF_E_TOPO_MISSING_SOURCE},
        {L"The topology loader does not support sink activates on output nodes.", MF_E_TOPO_SINK_ACTIVATES_UNSUPPORTED},
        {L"The topology loader does not support sink activates on output nodes.", MF_E_TOPO_SINK_ACTIVATES_UNSUPPORTED},
        {L"The topology loader does not support sink activates on output nodes.", MF_E_TOPO_SINK_ACTIVATES_UNSUPPORTED},
        {L"The sequencer cannot find a segment with the given ID.\n.", MF_E_SEQUENCER_UNKNOWN_SEGMENT_ID},
        {L"The context was canceled.\n.", MF_S_SEQUENCER_CONTEXT_CANCELED},
        {L"Cannot find source in source cache.\n.", MF_E_NO_SOURCE_IN_CACHE},
        {L"Cannot update topology flags.\n.", MF_S_SEQUENCER_SEGMENT_AT_END_OF_STREAM},
        {L"Cannot update topology flags.\n.", MF_S_SEQUENCER_SEGMENT_AT_END_OF_STREAM},
        {L"Cannot update topology flags.\n.", MF_S_SEQUENCER_SEGMENT_AT_END_OF_STREAM},
        {L"A valid type has not been set for this stream or a stream that it depends on.", MF_E_TRANSFORM_TYPE_NOT_SET},
        {L"A stream change has occurred. Output cannot be produced until the streams have been renegotiated.", MF_E_TRANSFORM_STREAM_CHANGE},
        {L"The transform cannot take the requested action until all of the input data it currently holds is processed or flushed.", MF_E_TRANSFORM_INPUT_REMAINING},
        {L"The transform requires a profile but no profile was supplied or found.", MF_E_TRANSFORM_PROFILE_MISSING},
        {L"The transform requires a profile but the supplied profile was invalid or corrupt.", MF_E_TRANSFORM_PROFILE_INVALID_OR_CORRUPT},
        {L"The transform requires a profile but the supplied profile ended unexpectedly while parsing.", MF_E_TRANSFORM_PROFILE_TRUNCATED},
        {L"The property ID does not match any property supported by the transform.", MF_E_TRANSFORM_PROPERTY_PID_NOT_RECOGNIZED},
        {L"The variant does not have the type expected for this property ID.", MF_E_TRANSFORM_PROPERTY_VARIANT_TYPE_WRONG},
        {L"An attempt was made to set the value on a read-only property.", MF_E_TRANSFORM_PROPERTY_NOT_WRITEABLE},
        {L"The array property value has an unexpected number of dimensions.", MF_E_TRANSFORM_PROPERTY_ARRAY_VALUE_WRONG_NUM_DIM},
        {L"The array or blob property value has an unexpected size.", MF_E_TRANSFORM_PROPERTY_VALUE_SIZE_WRONG},
        {L"The property value is out of range for this transform.", MF_E_TRANSFORM_PROPERTY_VALUE_OUT_OF_RANGE},
        {L"The property value is incompatible with some other property or mediatype set on the transform.", MF_E_TRANSFORM_PROPERTY_VALUE_INCOMPATIBLE},
        {L"The requested operation is not supported for the currently set output mediatype.", MF_E_TRANSFORM_NOT_POSSIBLE_FOR_CURRENT_OUTPUT_MEDIATYPE},
        {L"The requested operation is not supported for the currently set input mediatype.", MF_E_TRANSFORM_NOT_POSSIBLE_FOR_CURRENT_INPUT_MEDIATYPE},
        {L"The requested operation is not supported for the currently set combination of mediatypes.", MF_E_TRANSFORM_NOT_POSSIBLE_FOR_CURRENT_MEDIATYPE_COMBINATION},
        {L"The requested feature is not supported in combination with some other currently enabled feature.", MF_E_TRANSFORM_CONFLICTS_WITH_OTHER_CURRENTLY_ENABLED_FEATURES},
        {L"The transform cannot produce output until it gets more input samples.", MF_E_TRANSFORM_NEED_MORE_INPUT},
        {L"The requested operation is not supported for the current speaker configuration.", MF_E_TRANSFORM_NOT_POSSIBLE_FOR_CURRENT_SPKR_CONFIG},
        {L"The transform cannot accept mediatype changes in the middle of processing.", MF_E_TRANSFORM_CANNOT_CHANGE_MEDIATYPE_WHILE_PROCESSING},
        {L"The caller should not propagate this event to downstream components.", MF_S_TRANSFORM_DO_NOT_PROPAGATE_EVENT},
        {L"The input type is not supported for D3D device.", MF_E_UNSUPPORTED_D3D_TYPE},
        {L"The caller does not appear to support this transform's asynchronous capabilities.", MF_E_TRANSFORM_ASYNC_LOCKED},
        {L"An audio compression manager driver could not be initialized by the transform.", MF_E_TRANSFORM_CANNOT_INITIALIZE_ACM_DRIVER},
        {L"An audio compression manager driver could not be initialized by the transform.", MF_E_TRANSFORM_CANNOT_INITIALIZE_ACM_DRIVER},
        {L"An audio compression manager driver could not be initialized by the transform.", MF_E_TRANSFORM_CANNOT_INITIALIZE_ACM_DRIVER},
        {L"You are not allowed to open this file. Contact the content provider for further assistance.", MF_E_LICENSE_INCORRECT_RIGHTS},
        {L"The license for this media file has expired. Get a new license or contact the content provider for further assistance.", MF_E_LICENSE_OUTOFDATE},
        {L"You need a license to perform the requested operation on this media file.", MF_E_LICENSE_REQUIRED},
        {L"The licenses for your media files are corrupted. Contact Microsoft product support.", MF_E_DRM_HARDWARE_INCONSISTENT},
        {L"The APP needs to provide IMFContentProtectionManager callback to access the protected media file.", MF_E_NO_CONTENT_PROTECTION_MANAGER},
        {L"Client does not have rights to restore licenses.", MF_E_LICENSE_RESTORE_NO_RIGHTS},
        {L"Licenses are restricted and hence can not be backed up.", MF_E_BACKUP_RESTRICTED_LICENSE},
        {L"License restore requires machine to be individualized.", MF_E_LICENSE_RESTORE_NEEDS_INDIVIDUALIZATION},
        {L"Protection for stream is not required.", MF_S_PROTECTION_NOT_REQUIRED},
        {L"Component is revoked.", MF_E_COMPONENT_REVOKED},
        {L"Trusted functionality is currently disabled on this component.", MF_E_TRUST_DISABLED},
        {L"No Action is set on WMDRM Output Trust Authority.", MF_E_WMDRMOTA_NO_ACTION},
        {L"Action is already set on WMDRM Output Trust Authority.", MF_E_WMDRMOTA_ACTION_ALREADY_SET},
        {L"DRM Heaader is not available.", MF_E_WMDRMOTA_DRM_HEADER_NOT_AVAILABLE},
        {L"Current encryption scheme is not supported.", MF_E_WMDRMOTA_DRM_ENCRYPTION_SCHEME_NOT_SUPPORTED},
        {L"Action does not match with current configuration.", MF_E_WMDRMOTA_ACTION_MISMATCH},
        {L"Invalid policy for WMDRM Output Trust Authority.", MF_E_WMDRMOTA_INVALID_POLICY},
        {L"The policies that the Input Trust Authority requires to be enforced are unsupported by the outputs.", MF_E_POLICY_UNSUPPORTED},
        {L"The OPL that the license requires to be enforced are not supported by the Input Trust Authority.", MF_E_OPL_NOT_SUPPORTED},
        {L"The topology could not be successfully verified.", MF_E_TOPOLOGY_VERIFICATION_FAILED},
        {L"Signature verification could not be completed successfully for this component.", MF_E_SIGNATURE_VERIFICATION_FAILED},
        {L"Running this process under a debugger while using protected content is not allowed.", MF_E_DEBUGGING_NOT_ALLOWED},
        {L"MF component has expired.", MF_E_CODE_EXPIRED},
        {L"The current GRL on the machine does not meet the minimum version requirements.", MF_E_GRL_VERSION_TOO_LOW},
        {L"The current GRL on the machine does not contain any renewal entries for the specified revocation.", MF_E_GRL_RENEWAL_NOT_FOUND},
        {L"The current GRL on the machine does not contain any extensible entries for the specified extension GUID.", MF_E_GRL_EXTENSIBLE_ENTRY_NOT_FOUND},
        {L"The kernel isn't secure for high security level content.", MF_E_KERNEL_UNTRUSTED},
        {L"The response from protected environment driver isn't valid.", MF_E_PEAUTH_UNTRUSTED},
        {L"A non-PE process tried to talk to PEAuth.", MF_E_NON_PE_PROCESS},
        {L"We need to reboot the machine.", MF_E_REBOOT_REQUIRED},
        {L"Protection for this stream is not guaranteed to be enforced until the MEPolicySet event is fired.", MF_S_WAIT_FOR_POLICY_SET},
        {L"This video stream is disabled because it is being sent to an unknown software output.", MF_S_VIDEO_DISABLED_WITH_UNKNOWN_SOFTWARE_OUTPUT},
        {L"The GRL file is not correctly formed, it may have been corrupted or overwritten.", MF_E_GRL_INVALID_FORMAT},
        {L"The GRL file is in a format newer than those recognized by this GRL Reader.", MF_E_GRL_UNRECOGNIZED_FORMAT},
        {L"The GRL was reloaded and required all processes that can run protected media to restart.", MF_E_ALL_PROCESS_RESTART_REQUIRED},
        {L"The GRL was reloaded and the current process needs to restart.", MF_E_PROCESS_RESTART_REQUIRED},
        {L"The user space is untrusted for protected content play.", MF_E_USERMODE_UNTRUSTED},
        {L"PEAuth communication session hasn't been started.", MF_E_PEAUTH_SESSION_NOT_STARTED},
        {L"PEAuth's public key is revoked.", MF_E_PEAUTH_PUBLICKEY_REVOKED},
        {L"The GRL is absent.", MF_E_GRL_ABSENT},
        {L"The Protected Environment is trusted.", MF_S_PE_TRUSTED},
        {L"The Protected Environment is untrusted.", MF_E_PE_UNTRUSTED},
        {L"The Protected Environment Authorization service (PEAUTH) has not been started.", MF_E_PEAUTH_NOT_STARTED},
        {L"The sample protection algorithms supported by components are not compatible.", MF_E_INCOMPATIBLE_SAMPLE_PROTECTION},
        {L"No more protected environment sessions can be supported.", MF_E_PE_SESSIONS_MAXED},
        {L"WMDRM ITA does not allow protected content with high security level for this release.", MF_E_HIGH_SECURITY_LEVEL_CONTENT_NOT_ALLOWED},
        {L"WMDRM ITA cannot allow the requested action for the content as one or more components is not properly signed.", MF_E_TEST_SIGNED_COMPONENTS_NOT_ALLOWED},
        {L"WMDRM ITA does not support the requested action.", MF_E_ITA_UNSUPPORTED_ACTION},
        {L"WMDRM ITA encountered an error in parsing the Secure Audio Path parameters.", MF_E_ITA_ERROR_PARSING_SAP_PARAMETERS},
        {L"The Policy Manager action passed in is invalid.", MF_E_POLICY_MGR_ACTION_OUTOFBOUNDS},
        {L"The structure specifying Output Protection Level is not the correct format.", MF_E_BAD_OPL_STRUCTURE_FORMAT},
        {L"WMDRM ITA does not recognize the Explicite Analog Video Output Protection guid specified in the license.", MF_E_ITA_UNRECOGNIZED_ANALOG_VIDEO_PROTECTION_GUID},
        {L"IMFPMPHost object not available.", MF_E_NO_PMP_HOST},
        {L"WMDRM ITA could not initialize the Output Protection Level data.", MF_E_ITA_OPL_DATA_NOT_INITIALIZED},
        {L"WMDRM ITA does not recognize the Analog Video Output specified by the OTA.", MF_E_ITA_UNRECOGNIZED_ANALOG_VIDEO_OUTPUT},
        {L"WMDRM ITA does not recognize the Digital Video Output specified by the OTA.", MF_E_ITA_UNRECOGNIZED_DIGITAL_VIDEO_OUTPUT},
        {L"WMDRM ITA does not recognize the Digital Video Output specified by the OTA.", MF_E_ITA_UNRECOGNIZED_DIGITAL_VIDEO_OUTPUT},
        {L"WMDRM ITA does not recognize the Digital Video Output specified by the OTA.", MF_E_ITA_UNRECOGNIZED_DIGITAL_VIDEO_OUTPUT},
        {L"The continuity key supplied is not currently valid.", MF_E_CLOCK_INVALID_CONTINUITY_KEY},
        {L"No Presentation Time Source has been specified.", MF_E_CLOCK_NO_TIME_SOURCE},
        {L"The clock is already in the requested state.", MF_E_CLOCK_STATE_ALREADY_SET},
        {L"The clock has too many advanced features to carry out the request.", MF_E_CLOCK_NOT_SIMPLE},
        {L"Timer::SetTimer returns this success code if called happened while timer is stopped. Timer is not going to be dispatched until clock is running", MF_S_CLOCK_STOPPED},
        {L"Timer::SetTimer returns this success code if called happened while timer is stopped. Timer is not going to be dispatched until clock is running", MF_S_CLOCK_STOPPED},
        {L"Timer::SetTimer returns this success code if called happened while timer is stopped. Timer is not going to be dispatched until clock is running", MF_S_CLOCK_STOPPED},
        {L"The component does not support any more drop modes.", MF_E_NO_MORE_DROP_MODES},
        {L"The component does not support any more quality levels.", MF_E_NO_MORE_QUALITY_LEVELS},
        {L"The component does not support drop time functionality.", MF_E_DROPTIME_NOT_SUPPORTED},
        {L"Quality Manager needs to wait longer before bumping the Quality Level up.", MF_E_QUALITYKNOB_WAIT_LONGER},
        {L"Quality Manager is in an invalid state. Quality Management is off at this moment.", MF_E_QM_INVALIDSTATE},
        {L"Quality Manager is in an invalid state. Quality Management is off at this moment.", MF_E_QM_INVALIDSTATE},
        {L"Quality Manager is in an invalid state. Quality Management is off at this moment.", MF_E_QM_INVALIDSTATE},
        {L"No transcode output container type is specified.", MF_E_TRANSCODE_NO_CONTAINERTYPE},
        {L"The profile does not have a media type configuration for any selected source streams.", MF_E_TRANSCODE_PROFILE_NO_MATCHING_STREAMS},
        {L"Cannot find an encoder MFT that accepts the user preferred output type.", MF_E_TRANSCODE_NO_MATCHING_ENCODER},
        {L"Cannot find an encoder MFT that accepts the user preferred output type.", MF_E_TRANSCODE_NO_MATCHING_ENCODER},
        {L"Cannot find an encoder MFT that accepts the user preferred output type.", MF_E_TRANSCODE_NO_MATCHING_ENCODER},
        {L"Memory allocator is not initialized.", MF_E_ALLOCATOR_NOT_INITIALIZED},
        {L"Memory allocator is not committed yet.", MF_E_ALLOCATOR_NOT_COMMITED},
        {L"Memory allocator has already been committed.", MF_E_ALLOCATOR_ALREADY_COMMITED},
        {L"An error occurred in media stream.", MF_E_STREAM_ERROR},
        {L"Stream is not in a state to handle the request.", MF_E_INVALID_STREAM_STATE}
    };

    for (int i = 0; i < _countof(MFErrorMap); i++)
    {
        if (MFErrorMap[i].hr == hr)
        {
            errmsg.assign(MFErrorMap[i].wszMFError);
            return true;
        }
    }

    return false;
}


#if defined(VT_WINRT)
// we cannot use LoadLibrary in Win8, but can assume MF is available:
#pragma comment(lib, "mfplat")
#pragma comment(lib, "mfreadwrite")

// replacement for WinRT verboten MF API
// http://msdn.microsoft.com/en-us/library/windows/desktop/aa370467(v=vs.85).aspx
HRESULT STDAPICALLTYPE FrameRateToAverageTimePerFrame(UINT32 numerator, UINT32 denominator,
    UINT64* atpf)
{
    if (nullptr == atpf)
        return E_POINTER;

    *atpf = UINT64(.5 + 1e7 * double(denominator) / double(numerator));

    return S_OK;
}

// replacement for WinRT verboten MF API
// http://msdn.microsoft.com/en-us/library/windows/desktop/aa473720(v=vs.85).aspx
HRESULT STDAPICALLTYPE GetStrideForBitmapInfoHeader(DWORD format, DWORD dwWidth,
    LONG* pStride)
{
    if (nullptr == pStride)
        return E_POINTER;

#pragma message(__WARNING__"ignoring pixel format in stride computation")
    UNREFERENCED_PARAMETER(format);
    *pStride = LONG(4 * dwWidth);

    return S_OK;
}

#endif

//+----------------------------------------------------------------------------
//
// Class: CMediaFoundationAPI
// 
//-----------------------------------------------------------------------------
class CMediaFoundationAPI
{
public:
    CMediaFoundationAPI();

    bool IsAvailable() 
    { return pfnMFCreateSample != NULL; }

public:
    typedef HRESULT (STDAPICALLTYPE* T_MFStartup)(ULONG, DWORD); 
    typedef HRESULT (STDAPICALLTYPE* T_MFShutdown)();
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateAttributes)
        ( __out IMFAttributes**, __in UINT32);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateMediaType)
        (__deref_out IMFMediaType**);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateSourceResolver) 
        (__out IMFSourceResolver**);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateSourceReaderFromByteStream)
        (__in IMFByteStream*, __in_opt IMFAttributes*, __out IMFSourceReader**);
    typedef HRESULT (STDAPICALLTYPE* T_MFGetStrideForBitmapInfoHeader)
        (DWORD, DWORD, __out LONG*);
    typedef HRESULT (STDAPICALLTYPE* T_MFFrameRateToAverageTimePerFrame)
        (UINT32, UINT32, __out UINT64*);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateSinkWriterFromURL)
        (__in_opt LPCWSTR, __in_opt IMFByteStream*, __in_opt IMFAttributes*, __out IMFSinkWriter**);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateMemoryBuffer)
        (DWORD, __out IMFMediaBuffer**);
    typedef HRESULT (STDAPICALLTYPE* T_MFCreateSample)
        (__out IMFSample**);

public:
    T_MFStartup              pfnMFStartup;
    T_MFShutdown             pfnMFShutdown;
    T_MFCreateAttributes     pfnMFCreateAttributes;
    T_MFCreateMediaType      pfnMFCreateMediaType;
    T_MFCreateSourceResolver pfnMFCreateSourceResolver;
    T_MFCreateSourceReaderFromByteStream pfnMFCreateSourceReaderFromByteStream;
    T_MFGetStrideForBitmapInfoHeader     pfnMFGetStrideForBitmapInfoHeader;
    T_MFFrameRateToAverageTimePerFrame  pfnMFFrameRateToAverageTimePerFrame;
    T_MFCreateSinkWriterFromURL  pfnMFCreateSinkWriterFromURL;
    T_MFCreateMemoryBuffer  pfnMFCreateMemoryBuffer;
    T_MFCreateSample  pfnMFCreateSample;
};

CMediaFoundationAPI::CMediaFoundationAPI() :
   pfnMFStartup(NULL),
   pfnMFShutdown(NULL),
   pfnMFCreateAttributes(NULL),
   pfnMFCreateMediaType(NULL),
   pfnMFCreateSourceResolver(NULL),
   pfnMFCreateSourceReaderFromByteStream(NULL),
   pfnMFGetStrideForBitmapInfoHeader(NULL),
   pfnMFFrameRateToAverageTimePerFrame(NULL),
   pfnMFCreateSinkWriterFromURL(NULL),
   pfnMFCreateMemoryBuffer(NULL),
   pfnMFCreateSample(NULL)
{
#if defined(VT_WINRT)
    #define SET_MF_PROC(n) (pfn##n = ##n)
    
    SET_MF_PROC(MFStartup);
    SET_MF_PROC(MFShutdown);
    SET_MF_PROC(MFCreateAttributes);
    SET_MF_PROC(MFCreateMediaType);
    SET_MF_PROC(MFCreateSourceResolver);
    SET_MF_PROC(MFCreateSourceReaderFromByteStream);
    pfnMFGetStrideForBitmapInfoHeader = GetStrideForBitmapInfoHeader;
    pfnMFFrameRateToAverageTimePerFrame = FrameRateToAverageTimePerFrame;
    SET_MF_PROC(MFCreateSinkWriterFromURL);
    SET_MF_PROC(MFCreateMemoryBuffer);
    SET_MF_PROC(MFCreateSample);

    #undef SET_MF_PROC
#else
    HMODULE hpl = LoadLibraryW(L"mfplat.dll");
    if( !hpl )
    {    
        return;
    }

    HMODULE hrw = LoadLibraryW(L"mfreadwrite.dll");
    if( !hrw )
    {    
        return;
    }
    
    #define GET_MF_PROC(h, n)\
        if( (pfn##n = (T_##n)GetProcAddress(h, #n))==NULL ) return;
    
    GET_MF_PROC(hpl, MFStartup);
    GET_MF_PROC(hpl, MFShutdown);
    GET_MF_PROC(hpl, MFCreateAttributes);
    GET_MF_PROC(hpl, MFCreateMediaType);
    GET_MF_PROC(hpl, MFCreateSourceResolver);
    GET_MF_PROC(hrw, MFCreateSourceReaderFromByteStream);
    GET_MF_PROC(hpl, MFGetStrideForBitmapInfoHeader);
    GET_MF_PROC(hpl, MFFrameRateToAverageTimePerFrame);
    GET_MF_PROC(hrw, MFCreateSinkWriterFromURL);
    GET_MF_PROC(hpl, MFCreateMemoryBuffer);
    GET_MF_PROC(hpl, MFCreateSample);

    #undef GET_MF_PROC
#endif
};

CMediaFoundationAPI g_MediaFoundationAPI;

//#define DEBUG_METADATA
//#define DEBUG_LOG_ATTRIBUTES

//=============================================================================
// Global functions
//=============================================================================

//-----------------------------------------------------------------------------
// VtIsMediaFoundationAvailable
//-----------------------------------------------------------------------------

bool vt::VtIsMediaFoundationAvailable()
{ return g_MediaFoundationAPI.IsAvailable(); }

//-----------------------------------------------------------------------------
// VtCreateVideoSrc
//-----------------------------------------------------------------------------

HRESULT vt::VtCreateVideoSrc(IVideoSrc **ppVideoSrc, const WCHAR *filename) 
{
    HRESULT hr = S_OK;

    if (ppVideoSrc == NULL)
    {
        return E_POINTER;
    }
    
    if (!vt::VtIsMediaFoundationAvailable())
    {
        return E_NOINTERFACE;
    }

    IVideoSrc *videoSrc = VT_NOTHROWNEW CMFVideoSrc();
    if (!videoSrc)
    {
        return E_OUTOFMEMORY;
    }
    else
    {
        if (filename)
        {
            hr = videoSrc->OpenFile(filename);
            if (FAILED(hr))
            {
                videoSrc->Release();
                videoSrc = NULL;
            }
        }        
    }

    *ppVideoSrc = videoSrc;
    return hr;
}

#if defined(VT_WINRT)
HRESULT vt::VtCreateVideoSrc(IVideoSrc **ppVideoSrc, IStorageFile^ storageFile)
{
    HRESULT hr = S_OK;

    if (ppVideoSrc == NULL)
    {
        return E_POINTER;
    }
    
    if (!vt::VtIsMediaFoundationAvailable())
    {
        return E_NOINTERFACE;
    }

    IVideoSrc *videoSrc = VT_NOTHROWNEW CMFVideoSrc();
    if (!videoSrc)
    {
        return E_OUTOFMEMORY;
    }
    else
    {
        if (storageFile)
        {
            hr = videoSrc->OpenFile(storageFile);
            if (FAILED(hr))
            {
                videoSrc->Release();
                videoSrc = NULL;
            }
        }        
    }

    *ppVideoSrc = videoSrc;
    return hr;
}
#endif

//=============================================================================
// Constants
//=============================================================================

static const LONGLONG    MAX_FRAMES_TO_SKIP        = 100;
static const LONGLONG    TIME_HNS_TO_S_FACTOR    = 10000000;
/// static const GUID        FRAME_FORMAT            = MFVideoFormat_RGB32;
static const int        DEFAULT_IMAGE_OUT_TYPE    = OBJ_RGBIMG;
static const UCHAR        ALPHA_MAX_VALUE            = 255;

//=============================================================================
// Forward declarations
//=============================================================================

static HRESULT GetDefaultStride(IMFMediaType *pType, LONG *plStride);
static HRESULT SetAlphaChannel(vt::CRGBAImg &img, UCHAR value);
static HRESULT CreateMediaSource(PCWSTR sURL, IMFMediaSource **ppSource, 
                                 IMFByteStream **ppByteStream);

#ifdef DEBUG_METADATA
static HRESULT DumpMetadata(const WCHAR* wszFileName);
#endif // DEBUG_METADATA

#ifdef DEBUG_LOG_ATTRIBUTES
static HRESULT LogAllMediaTypes(IMFSourceReader* pReader);
static HRESULT LogAttributes(IMFAttributes *pAttributes);
static LPCWSTR GetGUIDNameConst(const GUID& guid);
static HRESULT GetGUIDName(const GUID& guid, WCHAR **ppwsz);
static HRESULT LogAttributeValueByIndex(IMFAttributes *pAttr, DWORD index);
static HRESULT SpecialCaseAttributeValue(GUID guid, const PROPVARIANT& var);
static void DBGMSG(PCWSTR format, ...);
#endif // DEBUG_LOG_ATTRIBUTES

//=============================================================================
// Inline functions
//=============================================================================

template <class T> void SafeRelease(T **ppT)
{
    if (*ppT)
    {
        (*ppT)->Release();
        *ppT = NULL;
    }
}

//=============================================================================
// CBufferLock class
//=============================================================================

class CBufferLock
{
public:
    CBufferLock(IMFMediaBuffer *pBuffer) : m_p2DBuffer(NULL), m_bLocked(FALSE)
    {
        m_pBuffer = pBuffer;
        m_pBuffer->AddRef();

        m_pBuffer->QueryInterface(IID_IMF2DBuffer, (void**)&m_p2DBuffer);
    }

    ~CBufferLock()
    {
        UnlockBuffer();
        SafeRelease(&m_pBuffer);
        SafeRelease(&m_p2DBuffer);
    }

    HRESULT LockBuffer(
        LONG  lDefaultStride,    // Minimum stride (with no padding).
        DWORD dwHeightInPixels,  // Height of the image, in pixels.
        BYTE  **ppbScanLine0,    // Receives a pointer to the first pixel in memory.
        LONG  *plStride          // Receives the actual stride.
        )
    {
        HRESULT hr = S_OK;

        // Use the 2-D version if available.
        if (m_p2DBuffer)
        {
            hr = m_p2DBuffer->Lock2D(ppbScanLine0, plStride);
            if (*plStride < 0)
            {
                // Lock2D always returns a pointer to the top row of pixels.
                // When stride is negative, the image is stored from bottom to top
                // and the top row is the last row in memory.  We want a pointer
                // to the first row in memory (the bottom row of pixels).
                *ppbScanLine0 += *plStride * (dwHeightInPixels - 1);
            }
        }
        else
        {
            // Use non-2D version.
            BYTE *pData = NULL;
            DWORD cbMaxLength = 0;
            DWORD cbCurrentLength = 0;
            hr = m_pBuffer->Lock(&pData, &cbMaxLength, &cbCurrentLength);
            if (SUCCEEDED(hr))
            {
                // Return the actual stride and a pointer to the start of the buffer.
                *plStride = lDefaultStride;
                *ppbScanLine0 = pData;
            }
        }

        m_bLocked = (SUCCEEDED(hr));

        return hr;
    }
    
    void UnlockBuffer()
    {
        if (m_bLocked)
        {
            if (m_p2DBuffer)
            {
                (void)m_p2DBuffer->Unlock2D();
            }
            else
            {
                (void)m_pBuffer->Unlock();
            }
            m_bLocked = FALSE;
        }
    }

private:
    IMFMediaBuffer  *m_pBuffer;
    IMF2DBuffer     *m_p2DBuffer;

    BOOL            m_bLocked;
};

//=============================================================================
// GetD3DTextureFromMFSample - 
//=============================================================================

#if defined(ENABLE_HARDWARE_DECODER)

HRESULT GetD3DTextureFromMFSample(IMFSample* pSample, 
    CComPtr<ID3D11Texture2D>& spTex2d)
{
    VT_HR_BEGIN()

    CComPtr<IMFMediaBuffer> spBuffer;
    VT_HR_EXIT( pSample->GetBufferByIndex(0, &spBuffer) );

    CComQIPtr<IMFDXGIBuffer> spDXGIBuffer(spBuffer);
    if (!spDXGIBuffer)
    {
        VT_HR_EXIT(E_NOINTERFACE);
    }

    UINT nSampleIndex = 0;
    VT_HR_EXIT( spDXGIBuffer->GetSubresourceIndex(&nSampleIndex) );
    //VT_ASSERT(nSampleIndex < MAX_SUBINDEX_TEXTURE);

    VT_HR_EXIT( spDXGIBuffer->GetResource(
            __uuidof(ID3D11Texture2D), 
            reinterpret_cast<void **>(&spTex2d)) );

    CComQIPtr<IDXGIResource> spDXGIResource(spTex2d);
    if (!spDXGIResource)
    {
        VT_HR_EXIT(E_NOINTERFACE);
    }

    // get descriptor
    D3D11_TEXTURE2D_DESC desc;
    spTex2d->GetDesc(&desc);

    VT_HR_END()
}

#endif

//=============================================================================
// CMFVideoSrc class
//=============================================================================

//-----------------------------------------------------------------------------
// CMFVideoSrc constructor
//-----------------------------------------------------------------------------

CMFVideoSrc::CMFVideoSrc() : m_bMFInitialized(false), m_pReader(NULL), 
    m_bCanSeek(false), m_iRefCount(1), m_eVideoFormat(VFNone),
    m_resizeWidth(-1), m_resizeHeight(-1), m_bSoftwareResize(false),
    m_rotationMultipleOf90(0)
{
    HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
    m_bComStarted = SUCCEEDED(hr);

    m_bTMStarted = false; 
	if( TaskManagerStartup() == S_OK )
	{
		m_bTMStarted = true;
	}
}


//-----------------------------------------------------------------------------
// CMFVideoSrc destructor
//-----------------------------------------------------------------------------

CMFVideoSrc::~CMFVideoSrc()
{
    Close();
    ShutdownMediaFoundation();
    if( m_bComStarted )
    {
        CoUninitialize();
    }
	if( m_bTMStarted )
	{
		TaskManagerShutdown();
	}
}

//-----------------------------------------------------------------------------
// OpenFile: Opens a video file.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::OpenFile(IN const WCHAR* wszFileName, 
                              VideoFormat eVideoFormat)
{
    IMFByteStream  *pByteStream = NULL;

    VT_HR_BEGIN()

    VT_HR_RET(InitializeMediaFoundation());

#ifdef DEBUG_METADATA
    DumpMetadata(wszFileName);
#endif // DEBUG_METADATA

    // Create a byte stream from the filename.
    VT_HR_EXIT( CreateMediaSource(wszFileName, NULL, &pByteStream) );

    // Open the video from the byte stream.
    VT_HR_EXIT( OpenByteStream(pByteStream, eVideoFormat) );

    // Store the current filename.
    VT_HR_EXIT( m_strCurrentFilename.assign(wszFileName) );

    VT_HR_EXIT_LABEL()
    SafeRelease(&pByteStream);
    
    return hr;
}

//-----------------------------------------------------------------------------
// OpenFile: Opens a video file.
//-----------------------------------------------------------------------------

#if defined(VT_WINRT)
HRESULT CMFVideoSrc::OpenFile(IStorageFile^ storageFile, VideoFormat eVideoFormat)
{
#if defined(MFVIDEOSRC_DECODE_TO_NV12)
    // always read NV12 from decoder for ARM devices, and convert to RGB as necessary
    eVideoFormat = NV12;
#endif

    IMFByteStream* pByteStream = NULL;

    VT_HR_BEGIN();

    try
    {
        VT_HR_EXIT(InitializeMediaFoundation());

        IRandomAccessStream^ stream = CSystem::Await(storageFile->OpenAsync(FileAccessMode::Read));

        // Create a byte stream from the random access stream.
        VT_HR_EXIT(MFCreateMFByteStreamOnStreamEx(reinterpret_cast<IUnknown *>(stream), &pByteStream));

        // Open the video from the byte stream.
        VT_HR_EXIT(OpenByteStream(pByteStream, eVideoFormat));
    }
    catch (...)
    {
        hr = E_FAIL;
    }

    VT_HR_EXIT_LABEL();
    SafeRelease(&pByteStream);

    return hr;
}
#endif

//-----------------------------------------------------------------------------
// OpenByteStream: Helper function used by both OpenFile overloads.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::OpenByteStream(IMFByteStream* pByteStream, VideoFormat eVideoFormat)
{
    IMFAttributes *pAttributes  = NULL;

    VT_HR_BEGIN()

    Close();

    m_eVideoFormat = eVideoFormat;
    switch(m_eVideoFormat)
    {
    case NV12:
        m_gFrameFormat = MFVideoFormat_NV12;
        break;
    case CIMG:
        m_gFrameFormat = MFVideoFormat_RGB32;
        break;
    case RGB32:
        m_gFrameFormat = MFVideoFormat_RGB32;
        break;
    default:
        return E_INVALIDARG;
        break;
    };

    // Configure the source reader to perform video processing.
    //
    // This includes:
    //   - YUV to RGB-32 (FRAME_FORMAT)
    //   - Software deinterlace
    VT_HR_EXIT( g_MediaFoundationAPI.pfnMFCreateAttributes(&pAttributes, 1 ) );
#if ( defined(ENABLE_HARDWARE_DECODER) && (_MSC_VER >= 1700) )
    // create Direct3D11 device and dxgi device manager to enable use of DXVA2
    D3D_FEATURE_LEVEL flret;
    m_pD3DDev = NULL;
    m_pD3DDevCon = NULL;
    hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 
        D3D11_CREATE_DEVICE_VIDEO_SUPPORT, 
        nullptr, 0, D3D11_SDK_VERSION, 
        &m_pD3DDev.p, &flret, &m_pD3DDevCon.p);

    if (hr == S_OK)
    {
        CComPtr<ID3D10Multithread> d3dmt;
        VT_HR_EXIT( m_pD3DDev->QueryInterface(IID_ID3D10Multithread, (void**)&d3dmt.p) );
        d3dmt->SetMultithreadProtected(TRUE);

        UINT resetToken;
        m_pdxgiDevMan = NULL;
        VT_HR_EXIT( MFCreateDXGIDeviceManager(&resetToken, &m_pdxgiDevMan.p) );
        VT_HR_EXIT( m_pdxgiDevMan->ResetDevice(m_pD3DDev, resetToken) );
        VT_HR_EXIT( pAttributes->SetUnknown(MF_SOURCE_READER_D3D_MANAGER, m_pdxgiDevMan) );

        VT_HR_EXIT( pAttributes->SetUINT32(MF_SOURCE_READER_DISABLE_DXVA, FALSE) );
#if (WINVER >= _WIN32_WINNT_WIN8)
        // this must be false when 'advanced video processing' is true
        VT_HR_EXIT( pAttributes->SetUINT32(MF_READWRITE_DISABLE_CONVERTERS, FALSE) ); 
        // enable GPU conversions
        VT_HR_EXIT( pAttributes->SetUINT32(MF_SOURCE_READER_ENABLE_ADVANCED_VIDEO_PROCESSING, TRUE) );
#endif
    } 
    else
#endif
    {
        VT_HR_EXIT( pAttributes->SetUINT32(MF_SOURCE_READER_ENABLE_VIDEO_PROCESSING, TRUE) );
        VT_HR_EXIT( pAttributes->SetUINT32(MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS, TRUE) );
    }

    // Create the source reader.
    VT_HR_EXIT( g_MediaFoundationAPI.pfnMFCreateSourceReaderFromByteStream(
        pByteStream, pAttributes, &m_pReader) );

#ifdef DEBUG_LOG_ATTRIBUTES
    LogAllMediaTypes(m_pReader);
#endif // DEBUG_LOG_ATTRIBUTES

    // Read presentation attributes.
    VT_HR_EXIT( ReadPresentationAttributes() );

    // Attempt to find a video stream.
    VT_HR_EXIT( SelectVideoStream() ); 

    // Read media type attributes.
    VT_HR_EXIT( ReadMediaTypeAttributes() );
    m_origWidth = m_videoFrameData.width;
    m_origHeight = m_videoFrameData.height;

    VT_HR_EXIT_LABEL()

    if (FAILED(hr))
    {
        Close();
    }
    SafeRelease(&pAttributes);
    
    return hr;
}

//---------------------------------------------------------------------------------------
// Close: Closes the current video file or stream.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::Close()
{
    HRESULT hr = S_OK;

    SafeRelease(&m_pReader);
    m_bCanSeek = false;
    m_videoFrameData = VIDEO_FRAME_DATA();
    m_strCurrentFilename.resize(0);
    return hr;
}


//-----------------------------------------------------------------------------
// Clone: 
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::Clone(IVideoSrc** ppClone)
{
    return vt::VtCreateVideoSrc(ppClone, m_strCurrentFilename.get_constbuffer());
}

//-----------------------------------------------------------------------------
// GetFileName: 
//-----------------------------------------------------------------------------
const WCHAR* CMFVideoSrc::GetFileName()
{
    return m_strCurrentFilename.size()? 
        m_strCurrentFilename.get_constbuffer(): NULL; 
}

//-----------------------------------------------------------------------------
// CanSeek: Indicates whether the current video file is seekable.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::CanSeek(OUT bool& canSeek)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    canSeek = m_bCanSeek;
    return S_OK;
}


//-----------------------------------------------------------------------------
// GetDuration: Provides the duration of the video file in seconds.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetPixelAspectRatio(double& pixelAspectRatio)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    pixelAspectRatio = m_videoFrameData.pixelAspectRatio;
    return S_OK;
}

//-----------------------------------------------------------------------------
// GetDuration: Provides the duration of the video file in seconds.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetDuration(double& durationInSeconds)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    durationInSeconds = m_videoFrameData.duration;
    return S_OK;
}

//-----------------------------------------------------------------------------
// GetFrameCount: Provides the number of frames in the video file.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetFrameCount(vt::Int32& frameCount) 
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    frameCount = vt::Int32(floor(m_videoFrameData.frameRate * m_videoFrameData.duration));
    return S_OK;
}

//-----------------------------------------------------------------------------
// GetFrameCount: Provides the frame rate of the video file.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetFrameRate(double& framesPerSecond)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    framesPerSecond = m_videoFrameData.frameRate;
    return S_OK;
}

HRESULT CMFVideoSrc::GetBitrate(float& bitsPerPixel)
{
    if (m_pReader == NULL || m_videoFrameData.width == 0 || m_videoFrameData.width == 0)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }
    bitsPerPixel = float(m_videoFrameData.bitrate) / float(m_videoFrameData.width) / float(m_videoFrameData.width);
    return S_OK;
}

HRESULT CMFVideoSrc::GetInterlaceMode(int& interlaceMode)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }
    interlaceMode = m_videoFrameData.interlaceMode;
    return S_OK;
}

//-----------------------------------------------------------------------------
// GetFrameSize: Provides the frame size of the video file.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetFrameSize(vt::Int32 &width, vt::Int32 &height)
{
    if (m_pReader == NULL) 
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    width  = m_videoFrameData.cropRect.Width();
    height = m_videoFrameData.cropRect.Height();

    // Adjust for non-square pixels.
    // Note that we change width (not height), even if a pixel aspect ratio
    // larger than one leads to upsampling.  This is to match the behavior
    // of Media Player and WPF's MediaElement.
    if (m_videoFrameData.pixelAspectRatio != 1)
    {
        width = vt::Int32(width * m_videoFrameData.pixelAspectRatio + 0.5);
    }

    return S_OK;
}

//-----------------------------------------------------------------------------
// GetFrameSize: Provides the frame size of the raw pixels in the video file.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::GetRawFrameSize(vt::Int32 &width, vt::Int32 &height)
{
    if (m_pReader == NULL) 
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    width  = m_videoFrameData.width;
    height = m_videoFrameData.height;

    return S_OK;
}

//-----------------------------------------------------------------------------
// Set* methods to control modifications to returned results
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::SetFrameSize(vt::Int32 width, vt::Int32 height, 
                                  bool bSoftwareResize)
{
    if (m_pReader == NULL) 
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }
    IMFMediaType *pType = NULL;
    VT_HR_BEGIN();

    VT_HR_EXIT( ((width>0)||(width==-1))?S_OK:E_INVALIDARG );
    if ((width > 0) && (!bSoftwareResize))
    {
        // width and height must be even
        VT_HR_EXIT( ((width&0x1)==1)?E_INVALIDARG:S_OK );
        VT_HR_EXIT( ((height&0x1)==1)?E_INVALIDARG:S_OK );
    }

#if defined(ENABLE_HARDWARE_RESIZE)
    bool bResetDXVAResize = ((width == -1) && (!m_bSoftwareResize));
#endif
    if (width == -1)
    {
        m_resizeWidth = -1;
        m_resizeHeight = -1;
        m_bSoftwareResize = false;
    }
    else
    {
        m_resizeWidth = width;
        m_resizeHeight = height;
#if defined(ENABLE_HARDWARE_RESIZE)
        m_bSoftwareResize = bSoftwareResize;
#else
        m_bSoftwareResize = true;
#endif
    }

#if defined(ENABLE_HARDWARE_RESIZE)
    if (!m_bSoftwareResize || bResetDXVAResize)
    {
        // reconfigure the source reader to set the frame size
        VT_HR_EXIT( g_MediaFoundationAPI.pfnMFCreateMediaType(&pType) );
        VT_HR_EXIT( pType->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video) );
        VT_HR_EXIT( pType->SetGUID(MF_MT_SUBTYPE, m_gFrameFormat) );
        if (bResetDXVAResize)
        {
            VT_HR_EXIT( MFSetAttributeSize(pType, MF_MT_FRAME_SIZE, m_origWidth, m_origHeight) );
        }
        else
        {
            VT_HR_EXIT( MFSetAttributeSize(pType, MF_MT_FRAME_SIZE, m_resizeWidth, m_resizeHeight) );
        }
        VT_HR_EXIT( m_pReader->SetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, NULL, pType) );
        ReadMediaTypeAttributes();
    }
#endif

VT_HR_EXIT_LABEL()
    SafeRelease(&pType);
    return hr;
}

HRESULT CMFVideoSrc::SetFrameRotation(int multipleOf90)
{
    m_rotationMultipleOf90 = multipleOf90;
    return S_OK;
}

//-----------------------------------------------------------------------------
// Seek: Seeks to the specified time in the video file.
//-----------------------------------------------------------------------------
// CAUTION: Does not check to see if the position was set correctly.

HRESULT CMFVideoSrc::Seek(IN double desiredFrameTimeInSeconds)
{
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    PROPVARIANT var;
    PropVariantInit(&var);
    var.vt = VT_I8;
    var.hVal.QuadPart = LONGLONG(desiredFrameTimeInSeconds * TIME_HNS_TO_S_FACTOR);
    return m_pReader->SetCurrentPosition(GUID_NULL, var);
}

//-----------------------------------------------------------------------------
// GetNextFrame methods
//-----------------------------------------------------------------------------

template <typename T> HRESULT CMFVideoSrc::GetNextFrameInternal(
    T** ppimg, double* pActualFrameTimeInSeconds, VideoFormat fmt)
{
    HRESULT            hr = S_OK;
    IMFSample        *pSample = NULL;
    IMFMediaBuffer  *pBuffer = NULL;
    DWORD            dwFlags = 0;
    LONGLONG        sampleTime;

#if !defined(MFVIDEOSRC_DECODE_TO_NV12)
    if (m_eVideoFormat != fmt)
    {
        return FIO_VIDEO_E_INVALID_VIDEO_FORMAT;
    }
#else
    fmt;
#endif

    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    VT_HR_EXIT( m_pReader->ReadSample( (DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM,
        0, NULL, &dwFlags, NULL, &pSample) );
    
    if (dwFlags & MF_SOURCE_READERF_ENDOFSTREAM)
    {
        hr = FIO_VIDEO_E_END_OF_STREAM;
        goto Exit;
    }

    if (dwFlags & MF_SOURCE_READERF_CURRENTMEDIATYPECHANGED)
    {
        // The media type changed. Get the attributes of the new media type.
        VT_HR_EXIT( ReadMediaTypeAttributes() );
    }

    if (pSample)
    {
        // Don't exit if it wasn't possible to get the time stamp of the current sample
        /*VT_HR_CHK*/ ( pSample->GetSampleTime(&sampleTime) );
        if (SUCCEEDED(hr) && pActualFrameTimeInSeconds != NULL)
        {
            *pActualFrameTimeInSeconds = double(sampleTime) / TIME_HNS_TO_S_FACTOR;
        }

#ifdef DEBUG_LOG_ATTRIBUTES
        DBGMSG(L"\n==> Sample at %g:\n", double(sampleTime) / TIME_HNS_TO_S_FACTOR);
        LogAttributes(pSample);
#endif // DEBUG_LOG_ATTRIBUTES

        // Create a contiguous buffer.
        VT_HR_RET( pSample->ConvertToContiguousBuffer(&pBuffer) );

        // Wrap the buffer in our own accessor class, then lock the buffer.
        CBufferLock lockedBuffer(pBuffer);
        BYTE *pBitmapData = NULL;
        LONG signedStride = 0;
        VT_HR_EXIT( lockedBuffer.LockBuffer(m_videoFrameData.stride, 
                                             m_videoFrameData.height, 
                                             &pBitmapData, &signedStride) );

        VT_HR_EXIT( this->ConvertSampleToImage(ppimg, pBitmapData,signedStride) );
    }
    else
    {
        hr = FIO_VIDEO_E_END_OF_STREAM;
    }

Exit:
    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
    return hr;
}

HRESULT CMFVideoSrc::GetNextFrame(IN OUT vt::CImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CImg* pimg = &img;
    return GetNextFrameInternal(&pimg, pActualFrameTimeInSeconds, CIMG);
}

HRESULT CMFVideoSrc::GetNextFrame(IN OUT vt::CRGB32VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CRGB32VideoImg* pimg = &img;
    return GetNextFrameInternal(&pimg, pActualFrameTimeInSeconds, RGB32);
}

HRESULT CMFVideoSrc::GetNextFrame(IN OUT vt::CNV12VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CNV12VideoImg* pimg = &img;
    return GetNextFrameInternal(&pimg, pActualFrameTimeInSeconds, NV12);
}

HRESULT CMFVideoSrc::GetNextFrame(IN OUT vt::IImageWriter* pDst, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDst, pDst };
    return GetNextFrameInternal(ppimgs, pActualFrameTimeInSeconds, m_eVideoFormat);
}

HRESULT CMFVideoSrc::GetNextFrame(IN OUT vt::IImageWriter* pDstY, IN OUT vt::IImageWriter* pDstUV, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDstY, pDstUV };
    return GetNextFrameInternal(ppimgs, pActualFrameTimeInSeconds, m_eVideoFormat);
}

//-----------------------------------------------------------------------------
// GetFrame methods
//-----------------------------------------------------------------------------

template <typename T> HRESULT CMFVideoSrc::GetFrameInternal(
    double desiredFrameTimeInSeconds, T** ppimg,
    double* pActualFrameTimeInSeconds, VideoFormat fmt)
{
    HRESULT     hr = S_OK;
    DWORD       dwFlags = 0;

    LONGLONG    hnsTimeStamp = 0;
    DWORD       cSkipped = 0;           // Number of skipped frames

    IMFSample *pSample = NULL;
    IMFMediaBuffer  *pBuffer = NULL;
    vt::Int64 hnsPos = vt::Int64(desiredFrameTimeInSeconds * TIME_HNS_TO_S_FACTOR);
    
#if !defined(MFVIDEOSRC_DECODE_TO_NV12)
    if (m_eVideoFormat != fmt)
    {
        return FIO_VIDEO_E_INVALID_VIDEO_FORMAT;
    }
#else
    fmt;
#endif
    if (m_pReader == NULL)
    {
        return FIO_VIDEO_E_NOT_INITIALIZED;
    }

    if (hnsPos < 0)
    {
        return E_INVALIDARG;
    }

    if (m_bCanSeek)
    {
        PROPVARIANT var;
        PropVariantInit(&var);

        var.vt = VT_I8;
        var.hVal.QuadPart = hnsPos;

        // Note: We subtract one tick because ReadSample will sometimes
        // return E_FAIL after seeking to the exact time of a frame in
        // the video.  This means we will nearly always seek to an
        // earlier frame in the video and read frames incrementally
        // until we get the desired frame.
        if (var.hVal.QuadPart > 0)
        {
            var.hVal.QuadPart--;
        }

        VT_HR_EXIT( m_pReader->SetCurrentPosition(GUID_NULL, var) );
    }
    else 
    {
        return FIO_VIDEO_E_BYTESTREAM_NOT_SEEKABLE;
    }

    // Pulls video frames from the source reader.

    // NOTE: Seeking might be inaccurate, depending on the container
    //       format and how the file was indexed. Therefore, the first
    //       frame that we get might be earlier than the desired time.
    //       If so, we skip up to MAX_FRAMES_TO_SKIP frames.

    while (true)
    {
        IMFSample *pSampleTmp = NULL;

        VT_HR_EXIT( m_pReader->ReadSample( (DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM,  0, NULL, &dwFlags, NULL, &pSampleTmp) );

        if (dwFlags & MF_SOURCE_READERF_ENDOFSTREAM)
        {
            break;
        }

        if (dwFlags & MF_SOURCE_READERF_CURRENTMEDIATYPECHANGED)
        {
            // The media type changed. Read the attributes of the new media type.
            VT_HR_EXIT( ReadMediaTypeAttributes() );
        }

        if (pSampleTmp == NULL)
        {
            continue;
        }

        // We got a sample. Hold onto it.

        SafeRelease(&pSample);

        pSample = pSampleTmp;
        pSample->AddRef();

        if ( SUCCEEDED( pSample->GetSampleTime(&hnsTimeStamp) ) )
        {
            // Keep going until we get a frame that is less than one frame from the
            // desired seek position.

            // During this process, we might reach the end of the file, so we
            // always cache the last sample that we got (pSample).
            if (hnsTimeStamp + m_videoFrameData.seekTolerance / 2 < hnsPos)
            {
                SafeRelease(&pSampleTmp);
                //wprintf(L"Asked for a sample at %ld ", hnsPos);
                //wprintf(L"Got %ld\n", hnsTimeStamp);
                ++cSkipped;
                continue;
            }
        }
        
        SafeRelease(&pSampleTmp);
        //wprintf(L"Number of skipped frames: %d\n", cSkipped);
        hnsPos = hnsTimeStamp;
        if (pActualFrameTimeInSeconds)
        {
            *pActualFrameTimeInSeconds = double(hnsPos) / TIME_HNS_TO_S_FACTOR;
        }
        break;
    }

    if (pSample)
    {
        // Create a contiguous buffer.
        VT_HR_RET( pSample->ConvertToContiguousBuffer(&pBuffer) );

        // Wrap the buffer in our own accessor class, then lock the buffer.
        CBufferLock lockedBuffer(pBuffer);
        BYTE *pBitmapData = NULL;
        LONG signedStride = 0;
        VT_HR_EXIT( lockedBuffer.LockBuffer(m_videoFrameData.stride, 
                                             m_videoFrameData.height, 
                                             &pBitmapData, &signedStride) );

        VT_HR_EXIT( this->ConvertSampleToImage(ppimg, pBitmapData,signedStride) );
    }
    else
    {
        hr = FIO_VIDEO_E_END_OF_STREAM;
    }

Exit:
    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
    return hr;
}

HRESULT CMFVideoSrc::GetFrame(IN double desiredFrameTimeInSeconds, 
    IN OUT vt::CImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CImg* pimg = &img;
    return GetFrameInternal(desiredFrameTimeInSeconds,&pimg,pActualFrameTimeInSeconds,CIMG);
}

HRESULT CMFVideoSrc::GetFrame(IN double desiredFrameTimeInSeconds, 
    IN OUT vt::CRGB32VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CRGB32VideoImg* pimg = &img;
    return GetFrameInternal(desiredFrameTimeInSeconds,&pimg,pActualFrameTimeInSeconds,RGB32);
}

HRESULT CMFVideoSrc::GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CNV12VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CNV12VideoImg* pimg = &img;
    return GetFrameInternal(desiredFrameTimeInSeconds,&pimg,pActualFrameTimeInSeconds,NV12);
}

HRESULT CMFVideoSrc::GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::IImageWriter* pDst, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDst, pDst };
    return GetFrameInternal(desiredFrameTimeInSeconds,ppimgs,pActualFrameTimeInSeconds,m_eVideoFormat);
}

HRESULT CMFVideoSrc::GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::IImageWriter* pDstY,  IN OUT vt::IImageWriter* pDstUV, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDstY, pDstUV };
    return GetFrameInternal(desiredFrameTimeInSeconds,ppimgs,pActualFrameTimeInSeconds,m_eVideoFormat);
}

template <typename T> HRESULT CMFVideoSrc::GetFrameInternal(
    Int32 frameIndex, T** ppimg,
    double* pActualFrameTimeInSeconds, VideoFormat fmt)
{
    HRESULT hr = S_OK;    
    double framesPerSecond = 0;
    double durationInSeconds = 0;

    VT_HR_EXIT( this->GetFrameRate(framesPerSecond) );
    VT_HR_EXIT( this->GetDuration(durationInSeconds) );
    vt::Int32 frameCount = vt::Int32(floor(framesPerSecond * durationInSeconds));

    if (frameIndex < 0 || frameIndex >= frameCount)
    {
        return E_INVALIDARG;
    }

    double desiredFrameTimeInSeconds = double(frameIndex) / framesPerSecond;
    VT_HR_EXIT( this->GetFrameInternal(desiredFrameTimeInSeconds, ppimg, pActualFrameTimeInSeconds, fmt) );

Exit:
    return hr;
}

HRESULT CMFVideoSrc::GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CImg* pimg = &img;
    return GetFrameInternal(frameIndex,&pimg,pActualFrameTimeInSeconds,CIMG);
}

HRESULT CMFVideoSrc::GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CRGB32VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CRGB32VideoImg* pimg = &img;
    return GetFrameInternal(frameIndex,&pimg,pActualFrameTimeInSeconds,RGB32);
}

HRESULT CMFVideoSrc::GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CNV12VideoImg& img, OUT double* pActualFrameTimeInSeconds)
{
    CNV12VideoImg* pimg = &img;
    return GetFrameInternal(frameIndex,&pimg,pActualFrameTimeInSeconds,NV12);
}

HRESULT CMFVideoSrc::GetFrame(IN vt::Int32 frameIndex, IN OUT vt::IImageWriter* pDst, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDst, pDst };
    return GetFrameInternal(frameIndex,ppimgs,pActualFrameTimeInSeconds,m_eVideoFormat);
}

HRESULT CMFVideoSrc::GetFrame(IN vt::Int32 frameIndex, IN OUT vt::IImageWriter* pDstY, IN OUT vt::IImageWriter* pDstUV, OUT double* pActualFrameTimeInSeconds)
{
    IImageWriter* ppimgs[2] = { pDstY, pDstUV };
    return GetFrameInternal(frameIndex,ppimgs,pActualFrameTimeInSeconds,m_eVideoFormat);
}

//-----------------------------------------------------------------------------
// AddRef
//-----------------------------------------------------------------------------

ULONG CMFVideoSrc::AddRef()
{
    return (ULONG)InterlockedIncrement(&m_iRefCount);
}

//-----------------------------------------------------------------------------
// Release
//-----------------------------------------------------------------------------

ULONG CMFVideoSrc::Release()
{
    if(InterlockedDecrement(&m_iRefCount)==0)
    {
        delete this;
        return 0;
    }
    return (ULONG)m_iRefCount;
}


//=============================================================================
// Private methods
//=============================================================================

//-----------------------------------------------------------------------------
// ReadPresentationAttributes: Processes attributes of the entire
// video (can seek, duration).
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::ReadPresentationAttributes()
{
    HRESULT hr = S_OK;
    PROPVARIANT var;
    PropVariantInit(&var);

    // Read media source characteristics to determine whether we can seek.
    // If the source has slow seeking, we will treat it as not supporting seeking.
    m_bCanSeek = false;
    VT_HR_EXIT( m_pReader->GetPresentationAttribute((DWORD)MF_SOURCE_READER_MEDIASOURCE,
        MF_SOURCE_READER_MEDIASOURCE_CHARACTERISTICS, &var) );
    VT_ASSERT( var.vt == VT_UI4 );
    UINT flags = var.uintVal;
    m_bCanSeek = (flags & MFMEDIASOURCE_CAN_SEEK) && !(flags & MFMEDIASOURCE_HAS_SLOW_SEEK);

    // Read the duration attribute.
    VT_HR_EXIT( m_pReader->GetPresentationAttribute((DWORD)MF_SOURCE_READER_MEDIASOURCE,
        MF_PD_DURATION, &var) );
    VT_ASSERT( var.vt == VT_UI8 );
    ULONGLONG durationInHundredsOfNanoseconds = var.uhVal.QuadPart;
    m_videoFrameData.duration = double(durationInHundredsOfNanoseconds) / TIME_HNS_TO_S_FACTOR;

Exit:
    // Clean up.
    (void)PropVariantClear(&var);
    return hr;
}

//-----------------------------------------------------------------------------
// SelectVideoStream: Finds the first video stream and sets the format.
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::SelectVideoStream()
{
    HRESULT hr = S_OK;
    IMFMediaType *pType = NULL;

    // Ensure only the first video stream is selected.
    VT_HR_EXIT( m_pReader->SetStreamSelection((DWORD)MF_SOURCE_READER_ALL_STREAMS, FALSE) );
    VT_HR_EXIT( m_pReader->SetStreamSelection((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, TRUE) );

    // Configure the source reader to return the created type.
    // The source reader will load the decoder if needed.
    VT_HR_EXIT( g_MediaFoundationAPI.pfnMFCreateMediaType(&pType) );
    VT_HR_EXIT( pType->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video) );
    VT_HR_EXIT( pType->SetGUID(MF_MT_SUBTYPE, m_gFrameFormat) );
    VT_HR_EXIT( m_pReader->SetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, NULL, pType) );
    
Exit:
    SafeRelease(&pType);
    return hr;
}

//-----------------------------------------------------------------------------

static MFOffset MakeOffset(float v)
{
    MFOffset offset;
    offset.value = short(v);
    offset.fract = WORD(65536 * (v - offset.value));
    return offset;
}

static MFVideoArea MakeArea(float x, float y, DWORD width, DWORD height)
{
    MFVideoArea area;
    area.OffsetX = MakeOffset(x);
    area.OffsetY = MakeOffset(y);
    area.Area.cx = width;
    area.Area.cy = height;
    return area;
}

static HRESULT GetVideoDisplayArea(IMFMediaType *pType, MFVideoArea *pArea, int &rectType)
{
    HRESULT hr = S_OK;

    // In pan-and-scan mode, try to get the pan-and-scan region.
    BOOL bPanScan = MFGetAttributeUINT32(pType, MF_MT_PAN_SCAN_ENABLED, FALSE);
    rectType = 2;
    if (bPanScan)
    {
        hr = pType->GetBlob(MF_MT_PAN_SCAN_APERTURE, (UINT8*)pArea, sizeof(MFVideoArea), NULL);
        rectType = 1;
    }

    // If not in pan-and-scan mode, or the pan-and-scan region is not set, 
    // get the minimimum display aperture.
    if (!bPanScan || hr == MF_E_ATTRIBUTENOTFOUND)
    {
        hr = pType->GetBlob(MF_MT_MINIMUM_DISPLAY_APERTURE, (UINT8*)pArea, sizeof(MFVideoArea), NULL);
        rectType = 2;

        if (hr == MF_E_ATTRIBUTENOTFOUND)
        {
            // Minimum display aperture is not set.

            // For backward compatibility with some components, 
            // check for a geometric aperture. 

            hr = pType->GetBlob(MF_MT_GEOMETRIC_APERTURE, (UINT8*)pArea, sizeof(MFVideoArea), NULL);
            if (SUCCEEDED(hr))
            {
                rectType = 3;
            }
        }

        // Default: Use the entire video area.

        if (hr == MF_E_ATTRIBUTENOTFOUND)
        {
            UINT32 width = 0, height = 0;
            hr = MFGetAttributeSize(pType, MF_MT_FRAME_SIZE, &width, &height);
            rectType = 2;

            if (SUCCEEDED(hr))
            {
                *pArea = MakeArea(0.0, 0.0, width, height);
            }
        }
    }
    return hr;
}

//-----------------------------------------------------------------------------
// ReadMediaTypeAttributes: Processes attributes of the current
// media type (frame rate, size, stride, and pixel aspect ratio).
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::ReadMediaTypeAttributes()
{
    HRESULT hr = S_OK;
    //GUID majorTypeGuid = GUID_NULL;
    GUID subtype = GUID_NULL;
    UINT32 numerator = 0;
    UINT32 denominator = 1;
    IMFMediaType *pType = NULL;

    // Get the first native media type from the stream.
    VT_HR_EXIT( m_pReader->GetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, &pType) );

#ifdef DEBUG_LOG_ATTRIBUTES
    DBGMSG(L"\n==> ReadMediaTypeAttributes:\n");
    LogAttributes(pType);
#endif // DEBUG_LOG_ATTRIBUTES

    // Make sure we got the expected video format.
    VT_HR_EXIT( pType->GetGUID(MF_MT_SUBTYPE, &subtype) );
    if (subtype != m_gFrameFormat)
    {
        VT_HR_EXIT(E_UNEXPECTED);
    }

    // Get the frame rate.
    hr = MFGetAttributeRatio(pType, MF_MT_FRAME_RATE, &numerator, &denominator);
    if ( (hr != S_OK) && (hr != MF_E_ATTRIBUTENOTFOUND) ) { VT_HR_EXIT(hr); }
    if (hr == S_OK)
    {
        m_videoFrameData.frameRate = double(numerator) / double(denominator);
        if (m_videoFrameData.frameRate > 0)
        {
            m_videoFrameData.seekTolerance = LONGLONG(TIME_HNS_TO_S_FACTOR/m_videoFrameData.frameRate);
        }
    }
    else
    {
        hr = S_OK;
    }

    // Get the frame size.
    VT_HR_EXIT( MFGetAttributeSize(pType, MF_MT_FRAME_SIZE, &m_videoFrameData.width, &m_videoFrameData.height) );
    // Check if DXVA resize is enabled and confirm that the resize is still
    // being applied - if not then switch to software resize
    // TODO: reapply DXVA resize?
    if ( (m_resizeWidth > 0) && (!m_bSoftwareResize) )
    {
        if ( ((int)m_videoFrameData.width != m_resizeWidth) ||
             ((int)m_videoFrameData.height != m_resizeHeight) )
        {
            m_bSoftwareResize = true;
        }
    }

    // Get the bitrate .
    if (FAILED(pType->GetUINT32(MF_MT_AVG_BITRATE, &m_videoFrameData.bitrate)) )
    {
        m_videoFrameData.bitrate = 0;
    }

    // Get the stride.
    VT_HR_EXIT( GetDefaultStride(pType, &m_videoFrameData.stride) );
    
    // Get the interlace mode
    VT_HR_EXIT( pType->GetUINT32(MF_MT_INTERLACE_MODE , &m_videoFrameData.interlaceMode) );

    // Get the display aperture.
    MFVideoArea videoArea;
    m_videoFrameData.cropRectType = 2;
    VT_HR_EXIT( GetVideoDisplayArea(pType, &videoArea, m_videoFrameData.cropRectType) );
    m_videoFrameData.cropRect.left = videoArea.OffsetX.value;
    m_videoFrameData.cropRect.top = videoArea.OffsetY.value;
    m_videoFrameData.cropRect.right = m_videoFrameData.cropRect.left + videoArea.Area.cx;
    m_videoFrameData.cropRect.bottom = m_videoFrameData.cropRect.top + videoArea.Area.cy;

    // Get the pixel aspect ratio.
    hr = MFGetAttributeRatio(pType, MF_MT_PIXEL_ASPECT_RATIO, &numerator, &denominator);
    if (SUCCEEDED(hr))
    {
        m_videoFrameData.pixelAspectRatio = double(numerator) / double(denominator);
    }
    else
    {
        // Assume square pixels if the pixel aspect ratio is unspecified.
        m_videoFrameData.pixelAspectRatio = 1;
        hr = S_OK;
    }

Exit:
    SafeRelease(&pType);
    return hr;
}

//-----------------------------------------------------------------------------
// ConvertSampleToImage: Converts a memory block of frame data to
// a Vision Tools image (CImg of DEFAULT_IMAGE_OUT_TYPE, C*VideoImg, or IImageWriter)
//
// pointer is always to start of image data; signedStride is negative if image
// needs to be flipped; the magnitude of signedStride is the stride
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::ConvertSampleToImage(vt::CImg** ppimg, BYTE* pBitmapData,
                                          LONG signedStride)
{
    VT_HR_BEGIN();
    // use 4 channel images for convert/wrap, resize, cropping, flipping, squarepix
    vt::CRGBAImg wrappedImage;
    vt::CRGBAImg resizedImage;
    vt::CRGBAImg croppedImage;
    vt::CRGBAImg flippedImage;
    vt::CRGBAImg sqarepixImage;

    int width = (int)m_videoFrameData.width;
    int height = (int)m_videoFrameData.height;
    LONG stride = abs(signedStride);

    if (m_eVideoFormat == NV12)
    {
        // decoder is returning NV12 image, so convert (sets alpha to 0xff)
        VT_HR_EXIT( wrappedImage.Create(width,height) );

        // wrap Y and UV data in CImg's for conversion
        vt::CLumaByteImg imgy;
        VT_HR_EXIT( imgy.Create(pBitmapData, width,height, stride) );
        vt::CUVByteImg imguv;
        VT_HR_EXIT( imguv.Create(pBitmapData+(height*stride), 
            width/2,height/2, stride) );

        // do convert
        VT_HR_EXIT( VtConvertImageNV12ToRGBA(wrappedImage, imgy, imguv) );
    }
    else
    {
        // decoder is returning RGB image, so just wrap in a CImg
        VT_HR_EXIT( wrappedImage.Create(pBitmapData, width, height, stride) );
        // Fill in the missing alpha channel - only when returning alpha.
        if ( ((**ppimg).Bands() == 4) ||
             (DEFAULT_IMAGE_OUT_TYPE == OBJ_RGBAIMG) )
        {
            VT_HR_EXIT( SetAlphaChannel(wrappedImage, ALPHA_MAX_VALUE) );
        }
    }
    const vt::CRGBAImg* pSourceImage = &wrappedImage;

    // Do sofware resize if enabled
    bool bResize = (m_bSoftwareResize && (m_resizeWidth > 0));
    if (bResize)
    {
        VT_HR_EXIT( resizedImage.Create(m_resizeWidth, m_resizeHeight) );
        VT_HR_EXIT( VtResizeImage(resizedImage, resizedImage.Rect(), 
            *pSourceImage, eSamplerKernelLanczos3) );
        pSourceImage = &resizedImage;
        width = m_resizeWidth;
        height = m_resizeHeight;
    }

    // If the frame has a display aperture, crop it.
    if (!bResize && (
            m_videoFrameData.cropRect.left  != 0 ||
            m_videoFrameData.cropRect.top   != 0 ||
            m_videoFrameData.cropRect.right  != width ||
            m_videoFrameData.cropRect.bottom != height) )
    {
        VT_HR_EXIT( pSourceImage->Share(croppedImage, &m_videoFrameData.cropRect) );
        pSourceImage = &croppedImage;
    }

    // Apply rotation and/or source flip
    if (m_rotationMultipleOf90&1) { std::swap(width,height); }
    if ( (m_rotationMultipleOf90 != 0) || (signedStride < 0) )
    {
        VT_HR_EXIT( flippedImage.Create(width,height) );
        VT_HR_EXIT( VtRotateImage(flippedImage, *pSourceImage, 
            m_rotationMultipleOf90+((signedStride<0)?(2):(0))) );
        pSourceImage = &flippedImage;
    }

    // Adjust for non-square pixels if needed.
    if (!bResize)
    {
        int srcw, srch;
        VT_HR_EXIT( GetFrameSize(srcw, srch) );
        if (m_rotationMultipleOf90&1) { std::swap(srcw,srch); }
        if (pSourceImage->Width() != srcw || pSourceImage->Height() != srch)
        {
            VT_HR_EXIT( vt::VtResizeImage(sqarepixImage, vt::CRect(0, 0, srcw, srch),
                                           *pSourceImage, eSamplerKernelLanczos3) );
            pSourceImage = &sqarepixImage;
        }
    }

    // Make sure the destination image is the right size, then convert.
    VT_HR_EXIT( vt::CreateImageForTransform(**ppimg, width, height,
                                            DEFAULT_IMAGE_OUT_TYPE) );
    VT_HR_EXIT( vt::VtConvertImage(**ppimg, *pSourceImage) );

    VT_HR_END();
}

CVideoImgInfo& CMFVideoSrc::VideoImgInfo(int w, int h)
{
    if ( w == -1)
    {
        m_VideoImgInfo.rectValidPixels = m_videoFrameData.cropRect;
    }
    else
    {
        m_VideoImgInfo.rectValidPixels = CRect(0,0,w,h);
    }
    // Debugging...
    // vi.rectValidPixels = CRect(20,20,m_videoFrameData.cropRect.Width()-20, m_videoFrameData.cropRect.Height()-20);
    switch(m_videoFrameData.cropRectType)
    {
    case 1:
        m_VideoImgInfo.eValidPixelsRectType = CVideoImgInfo::PanScan;
        break;
    case 2:
    default:
        m_VideoImgInfo.eValidPixelsRectType = CVideoImgInfo::MinDisplay;
        break;
    case 3:
        m_VideoImgInfo.eValidPixelsRectType = CVideoImgInfo::Geometric;
        break;
    };
    m_VideoImgInfo.dAspectRatio   = m_videoFrameData.pixelAspectRatio;
    m_VideoImgInfo.iInterlaceMode = m_videoFrameData.interlaceMode;
    return m_VideoImgInfo;
}

HRESULT CMFVideoSrc::ConvertSampleToImage(vt::CRGB32VideoImg** ppimg, BYTE* pBitmapData,
                                          LONG signedStride)
{
    VT_HR_BEGIN()

    // create dest image
    bool bResize = (m_bSoftwareResize && (m_resizeWidth > 0));
    int dstW = (bResize)?(m_resizeWidth):(m_videoFrameData.width);
    int dstH = (bResize)?(m_resizeHeight):(m_videoFrameData.height);
    if (m_rotationMultipleOf90&1) { std::swap(dstW,dstH); }
    VT_HR_EXIT( (**ppimg).Create(dstW, dstH, VideoImgInfo(dstW,dstH)) );

    // create readerwriter for dest image
    CImgReaderWriter<CRGBAByteImg> imgrw;
    VT_HR_EXIT( imgrw.Create(dstW,dstH,dstW,dstH,CPoint(0,0)) );
    VT_HR_EXIT( (**ppimg).GetImg().Share(imgrw) );

    VT_HR_EXIT( ProcessVideoFrame(&imgrw, nullptr, pBitmapData, signedStride) );

    VT_HR_END();
}

HRESULT CMFVideoSrc::ConvertSampleToImage(vt::CNV12VideoImg** ppimg, BYTE* pBitmapData,
                                          LONG signedStride)
{
    VT_HR_BEGIN()

    bool bResize = (m_bSoftwareResize && (m_resizeWidth > 0));
    int dstW = (bResize)?(m_resizeWidth):(m_videoFrameData.width);
    int dstH = (bResize)?(m_resizeHeight):(m_videoFrameData.height);
    if (m_rotationMultipleOf90&1) { std::swap(dstW,dstH); }
    VT_HR_EXIT( (**ppimg).Create(dstW, dstH, VideoImgInfo(dstW,dstH)) );

    // create readerwriter for dest  img's
    CImgReaderWriter<CLumaByteImg> imgrwY;
    VT_HR_EXIT( imgrwY.Create(dstW,dstH,dstW,dstH,CPoint(0,0)) );
    VT_HR_EXIT( (**ppimg).GetYImg().Share(imgrwY) );
    CImgReaderWriter<CUVByteImg> imgrwUV;
    VT_HR_EXIT( imgrwUV.Create(dstW/2,dstH/2,dstW/2,dstH/2,CPoint(0,0)) );
    VT_HR_EXIT( (**ppimg).GetUVImg().Share(imgrwUV) );

    VT_HR_EXIT( ProcessVideoFrame(&imgrwY, &imgrwUV, pBitmapData, signedStride) );

    VT_HR_END();
}

HRESULT CMFVideoSrc::ConvertSampleToImage(vt::IImageWriter** ppimgs, BYTE* pBitmapData,
                                          LONG signedStride)
{
    VT_HR_BEGIN();
    VT_HR_EXIT( ProcessVideoFrame(ppimgs[0], ppimgs[1], pBitmapData, signedStride) );
    VT_HR_END();
}

//=============================================================================
// Functions to build and run a transform graph for resize, rotation,
// NV12->RGB conversion, or just copy if none of these is needed
//=============================================================================

HRESULT CMFVideoSrc::SetProcessVideoFrameTransform(
    IImageWriter* pDst, IImageWriter* pDstUV,
    IImageReader* pSrc, IImageReader* pSrcUV,
    bool bFlip)
{
    VT_HR_BEGIN();

    bool bNV12Dst = (pDstUV != nullptr) && (pDstUV != pDst);
    bool bNV12Src = (pSrcUV != nullptr);

    int srcW = m_videoFrameData.width;
    int srcH = m_videoFrameData.height;

    int totalRot90 = (m_rotationMultipleOf90&3) + (bFlip?2:0);

    bool bResize = (m_resizeWidth > 0) && ( (m_resizeWidth != srcW) || (m_resizeHeight != srcH) );
    int dstW = bResize?m_resizeWidth:srcW;
    int dstH = bResize?m_resizeHeight:srcH;

    int nidx = 0; // fill node arrays as needed

    // set up rotation transforms and nodes
    if (totalRot90 != 0)
    {
        m_vpgraph.node[nidx].SetTransform(&m_vpgraph.xformRot[0]);
        if (bNV12Src)
        {
            m_vpgraph.xformRot[0].InitializeRotate(dstW,dstH,totalRot90,VT_IMG_FIXED(OBJ_LUMAIMG));
            m_vpgraph.xformRot[1].InitializeRotate(dstW/2,dstH/2,totalRot90,VT_IMG_FIXED(OBJ_UVIMG));
            m_vpgraph.nodeuv[nidx].SetTransform(&m_vpgraph.xformRot[1]);
        }
        else
        {
            m_vpgraph.xformRot[0].InitializeRotate(dstW,dstH,totalRot90,VT_IMG_FIXED(OBJ_RGBAIMG));
        }
        nidx++;
    }

    // set up resize transforms and nodes
    if (bResize)
    {
        float sx = (float)srcW/(float)m_resizeWidth;
        float sy = (float)srcH/(float)m_resizeHeight;
        m_vpgraph.node[nidx].SetTransform(&m_vpgraph.xformRsz[0]);
        if (bNV12Src)
        {
            VT_HR_EXIT( m_vpgraph.xformRsz[0].InitializeResize(sx,0.f,sy,0.f,VT_IMG_FIXED(OBJ_LUMAIMG),eSamplerKernelLanczos3) );
            VT_HR_EXIT( m_vpgraph.xformRsz[1].InitializeResize(sx,0.f,sy,0.f,VT_IMG_FIXED(OBJ_UVIMG),eSamplerKernelLanczos3) );
            m_vpgraph.nodeuv[nidx].SetTransform(&m_vpgraph.xformRsz[1]);
        }
        else
        {
            VT_HR_EXIT( m_vpgraph.xformRsz[0].InitializeResize(sx,0.f,sy,0.f,VT_IMG_FIXED(OBJ_RGBAIMG),eSamplerKernelLanczos3) );
        }
        nidx++;
    }

    // set pointers to top nodes and wire up graph
    //
    // can be:   ysrc[->Resize][->Rotate]->Convert->rgbdst
    //          uvsrc[->Resize][->Rotate]--^
    //     or:   ysrc[->Resize][->Rotate]->ydst
    //          uvsrc[->Resize][->Rotate]->uvdst
    //     or:   ysrc->Copy->ydst
    //          uvsrc->Copy->uvdst
    //     or: rgbsrc[->Resize][->Rotate]->rgbdst
    //     or: rgbsrc->Copy->rgbdst
    //
    if (!bNV12Dst) { m_vpgraph.pTopNodeUV = nullptr; }
    if (!bNV12Dst && bNV12Src)
    {
        // src is NV12 and dest is RGB so top node is convert
        VT_HR_EXIT( m_vpgraph.xformCvt.Initialize() );
        m_vpgraph.nodeCvt.SetTransform(&m_vpgraph.xformCvt);
        VT_HR_EXIT( m_vpgraph.nodeCvt.SetSourceCount(2) );
        m_vpgraph.pTopNode = &m_vpgraph.nodeCvt;

        // bind either sources or next nodes to convert node
        if (nidx > 0)
        {
            // rotate and/or resize are present, so bind the next nodes to convert
            VT_HR_EXIT( m_vpgraph.nodeCvt.BindSourceToTransform(0, &m_vpgraph.node[0]) );
            VT_HR_EXIT( m_vpgraph.nodeCvt.BindSourceToTransform(1, &m_vpgraph.nodeuv[0]) );
        }
        else
        {
            // no rotate or resize so bind sources to convert
            VT_HR_EXIT( m_vpgraph.nodeCvt.BindSourceToReader(0, pSrc, &m_vpgraph.ex) );
            VT_HR_EXIT( m_vpgraph.nodeCvt.BindSourceToReader(1, pSrcUV, &m_vpgraph.ex) );
        }
    }
    else
    {
        // no convert, so node[0(,1)] are top nodes, and can be copy, rotate, or resize
        if (nidx == 0)
        {
            // no rotation or resize, so initialize copy transform
            VT_HR_EXIT( m_vpgraph.xformCpy.Initialize(true) );
            m_vpgraph.node[0].SetTransform(&m_vpgraph.xformCpy);
            if (bNV12Dst) { m_vpgraph.nodeuv[0].SetTransform(&m_vpgraph.xformCpy); }
            nidx++;
        }
        m_vpgraph.pTopNode = &m_vpgraph.node[0];
        if (bNV12Dst) { m_vpgraph.pTopNodeUV = &m_vpgraph.nodeuv[0]; }
    }

    // bind source(s), if not already attached to convert, to last node(s)
    if (nidx > 0)
    {
        VT_HR_EXIT( m_vpgraph.node[nidx-1].BindSourceToReader(pSrc, &m_vpgraph.ex) );
        if (bNV12Src) { VT_HR_EXIT( m_vpgraph.nodeuv[nidx-1].BindSourceToReader(pSrcUV, &m_vpgraph.ex) ); }
    }

    // bind second level of nodes to first level if active
    if (nidx == 2)
    {
        VT_HR_EXIT( m_vpgraph.node[0].BindSourceToTransform(&m_vpgraph.node[1]) );
        if (bNV12Src) { VT_HR_EXIT( m_vpgraph.nodeuv[0].BindSourceToTransform(&m_vpgraph.nodeuv[1]) ); }
    }

    VT_HR_END();
}

HRESULT CMFVideoSrc::ProcessVideoFrame(IImageWriter* pDst, IImageWriter* pDstUV,
    BYTE* pSrc, int signedStride)
{
    VT_HR_BEGIN();

    int stride = abs(signedStride);
    int srcW = m_videoFrameData.width;
    int srcH = m_videoFrameData.height;

    // set up image source reader/writer and set up current transform graph
    CImgReaderWriter<CLumaByteImg> imgrwY;
    CImgReaderWriter<CUVByteImg> imgrwUV;
    CImgReaderWriter<CRGBAByteImg> imgrw;
    if (m_eVideoFormat == NV12)
    {
        VT_HR_EXIT( ((CLumaByteImg&)imgrwY).Create(pSrc,srcW,srcH,stride) );
        VT_HR_EXIT( ((CUVByteImg&)imgrwUV).Create(pSrc+(stride*srcH),srcW/2,srcH/2,stride) );
        SetProcessVideoFrameTransform(pDst,pDstUV,&imgrwY,&imgrwUV,signedStride<0);
    }
    else
    {
        VT_HR_EXIT( ((CRGBAByteImg&)imgrw).Create(pSrc,srcW,srcH,stride) );
        SetProcessVideoFrameTransform(pDst,pDstUV,&imgrw,nullptr,signedStride<0);
    }

    // invoke RGBA/Y graph
    {
        CRasterBlockMap bm(BLOCKITER_INIT(pDst->GetImgInfo().Rect(), 256, 256));
        m_vpgraph.pTopNode->SetDest(NODE_DEST_PARAMS(pDst, &bm));
        VT_HR_EXIT( PushTransformTaskAndWait(m_vpgraph.pTopNode, (CTaskStatus*)NULL) ); 
    }
    // invoke UV graph if necessary
    if ( m_vpgraph.pTopNodeUV != nullptr )
    {
        CRasterBlockMap bm(BLOCKITER_INIT(pDstUV->GetImgInfo().Rect(), 256, 256));
        m_vpgraph.pTopNodeUV->SetDest(NODE_DEST_PARAMS(pDstUV, &bm));
        VT_HR_EXIT( PushTransformTaskAndWait(m_vpgraph.pTopNodeUV, (CTaskStatus*)NULL) ); 
    }

    VT_HR_END();
}

//-----------------------------------------------------------------------------
// InitializeMediaFoundation
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::InitializeMediaFoundation()
{
    HRESULT hr = S_OK;
    if (!m_bMFInitialized)
    {
        if (!VtIsMediaFoundationAvailable())
        {
            return E_NOINTERFACE;
        }
        hr = g_MediaFoundationAPI.pfnMFStartup(MF_VERSION, MFSTARTUP_FULL);
        m_bMFInitialized = true;
    }
    return hr;
}

//-----------------------------------------------------------------------------
// ShutdownMediaFoundation
//-----------------------------------------------------------------------------

HRESULT CMFVideoSrc::ShutdownMediaFoundation()
{    
    HRESULT hr = S_OK;
    if (m_bMFInitialized)
    {    
        hr = g_MediaFoundationAPI.pfnMFShutdown();
        m_bMFInitialized = false;
    }
    return hr;
}


//=============================================================================
// Static functions
//=============================================================================

//-----------------------------------------------------------------------------
// GetDefaultStride: Gets the default stride for a given media type.
//-----------------------------------------------------------------------------

static HRESULT GetDefaultStride(IMFMediaType *pType, LONG *plStride)
{
    LONG lStride = 0;

    // Try to get the default stride from the media type.
    HRESULT hr = pType->GetUINT32(MF_MT_DEFAULT_STRIDE, (UINT32*)&lStride);
    if (FAILED(hr))
    {
        // Attribute not set. Try to calculate the default stride.

        // Get the subtype and the image size.
        GUID subtype = GUID_NULL;
        VT_HR_EXIT( pType->GetGUID(MF_MT_SUBTYPE, &subtype) );

        UINT32 width = 0;
        UINT32 height = 0;
        VT_HR_EXIT( MFGetAttributeSize(pType, MF_MT_FRAME_SIZE, &width, &height) );

        VT_HR_EXIT( g_MediaFoundationAPI.pfnMFGetStrideForBitmapInfoHeader(
            subtype.Data1, width, &lStride) );

        // Set the attribute for later reference.
        VT_HR_EXIT( pType->SetUINT32(MF_MT_DEFAULT_STRIDE, UINT32(lStride)) );
    }

    if (SUCCEEDED(hr))
    {
        *plStride = lStride;
    }

Exit:
    return hr;
}

//-----------------------------------------------------------------------------
// SetAlphaChannel: Fills the alpha channel of the given image.
//-----------------------------------------------------------------------------

static HRESULT SetAlphaChannel(vt::CRGBAImg &img, UCHAR value)
{
    // Note: Assumes that the image is valid.
    HRESULT hr = S_OK;
    INT i = 0, j = 0;

    for (i = 0; i < img.Height(); i++)
    {
        vt::RGBAPix* pPix = img.Ptr(i);
        for (j = 0; j < img.Width(); j++, pPix++)
        {
            pPix->a = value;
        }
    }

    return hr;
}

//-----------------------------------------------------------------------------
// CreateMediaSource: Creates a media source or a byte stream from
// the given URL.
//-----------------------------------------------------------------------------

static HRESULT CreateMediaSource(PCWSTR sURL, IMFMediaSource **ppSource, 
                                 IMFByteStream **ppByteStream)
{
    MF_OBJECT_TYPE objectType = MF_OBJECT_INVALID;

    IMFSourceResolver* pSourceResolver = NULL;
    IUnknown* pSource = NULL;
    DWORD createObjFlags = MF_RESOLUTION_MEDIASOURCE;
    // Create the source resolver.
    HRESULT hr = g_MediaFoundationAPI.pfnMFCreateSourceResolver(&pSourceResolver);
    if (FAILED(hr))
    {
        goto done;
    }

    // Use the source resolver to create the media source.
    if (ppSource == NULL)
    {
        createObjFlags = MF_RESOLUTION_BYTESTREAM;
    }

    createObjFlags |= MF_RESOLUTION_CONTENT_DOES_NOT_HAVE_TO_MATCH_EXTENSION_OR_MIME_TYPE;

    // Note: For simplicity this sample uses the synchronous method to create 
    // the media source. However, creating a media source can take a noticeable
    // amount of time, especially for a network source. For a more responsive 
    // UI, use the asynchronous BeginCreateObjectFromURL method.

    hr = pSourceResolver->CreateObjectFromURL(
        sURL,                       // URL of the source.
        createObjFlags,                // Create a source object.
        NULL,                       // Optional property store.
        &objectType,                // Receives the created object type. 
        &pSource                    // Receives a pointer to the media source.
        );
    if (FAILED(hr))
    {
        goto done;
    }

    // Get the IMFMediaSource interface from the media source.
    if (objectType == MF_OBJECT_MEDIASOURCE && ppSource != NULL)
    {
        hr = pSource->QueryInterface(IID_PPV_ARGS(ppSource));
    }
    else if (objectType == MF_OBJECT_BYTESTREAM && ppByteStream != NULL)
    {
        hr = pSource->QueryInterface(IID_PPV_ARGS(ppByteStream));
    }
    else
    {
        hr = E_UNEXPECTED;
    }

done:
    SafeRelease(&pSourceResolver);
    SafeRelease(&pSource);
    return hr;
}


//=============================================================================
// Debugging functions
//=============================================================================

#ifdef DEBUG_METADATA

static HRESULT DisplayMetadata(IMFPresentationDescriptor* pPD, 
                               IMFMetadataProvider* pProvider, DWORD dwStreamId)
{
    HRESULT hr = S_OK;
    IMFMetadata* pMetadata = NULL;
    PROPVARIANT varNames;

    VT_HR_EXIT( pProvider->GetMFMetadata(pPD, dwStreamId, 0, &pMetadata) );
    VT_HR_EXIT( pMetadata->GetAllPropertyNames(&varNames) );

    for (ULONG i = 0; i < varNames.calpwstr.cElems; i++)
    {
        DebuggerPrint("%S\n", varNames.calpwstr.pElems[i]);

        PROPVARIANT varValue;
        hr = pMetadata->GetProperty( varNames.calpwstr.pElems[i], &varValue );
        if (SUCCEEDED(hr))
        {
            //DisplayProperty(varValue);
            PropVariantClear(&varValue);
        }
    }

    PropVariantClear(&varNames);

Exit:
    SafeRelease(&pMetadata);
    return hr;
}

static HRESULT DumpMediaMetadata(IMFMediaSource* pSource)
{
    HRESULT hr = S_OK;
    IMFPresentationDescriptor* pPD = NULL;
    IMFGetService* pGetService = NULL;
    IMFMetadataProvider* pProvider = NULL;
    IMFStreamDescriptor* pDescriptor = NULL;

    hr = pSource->CreatePresentationDescriptor(&pPD);
    if (FAILED(hr))
    {
        DebuggerPrint("Media source can't create presentation descriptor.\n");
        goto Exit;
    }
    hr = pSource->QueryInterface(&pGetService);
    if (FAILED(hr))
    {
        DebuggerPrint("Media source doesn't support GetService.\n");
        goto Exit;
    }
    hr = pGetService->GetService(MF_METADATA_PROVIDER_SERVICE, IID_PPV_ARGS(&pProvider));
    if (FAILED(hr))
    {
        DebuggerPrint("Media source doesn't support metadata provider service.\n");
        goto Exit;
    }

    DebuggerPrint("--- Presentation metadata ---\n");
    DisplayMetadata(pPD, pProvider, 0);

    DWORD descriptorCount;
    VT_HR_EXIT( pPD->GetStreamDescriptorCount(&descriptorCount) );

    for (DWORD i = 0; i < descriptorCount; i++)
    {
        BOOL fIsSelected;
        VT_HR_EXIT( pPD->GetStreamDescriptorByIndex(i, &fIsSelected, &pDescriptor) );
        DWORD dwStreamId;
        VT_HR_EXIT( pDescriptor->GetStreamIdentifier(&dwStreamId) );
        SafeRelease(&pDescriptor);
        DebuggerPrint("--- Stream %d metadata (id = %d) ---\n", i, dwStreamId);
        DisplayMetadata(pPD, pProvider, 0);
    }

Exit:
    SafeRelease(&pPD);
    SafeRelease(&pGetService);
    SafeRelease(&pProvider);
    return hr;
}

static HRESULT DisplayProperty(PROPERTYKEY& key, PROPVARIANT& propvar)
{
    PWSTR pszValue;

    HRESULT hr = PSFormatForDisplayAlloc(key, propvar, PDFF_PREFIXNAME, &pszValue);
    if (SUCCEEDED(hr))
    {
        DebuggerPrint("%d\t%S\n", key.pid, pszValue);
        CoTaskMemFree(pszValue);
    }
    else
    {
        hr = PSFormatForDisplayAlloc(key, propvar, PDFF_DEFAULT, &pszValue);
        if (SUCCEEDED(hr))
        {
            DebuggerPrint("%d\tUnknown: %S\n", key.pid, pszValue);
            CoTaskMemFree(pszValue);
        }
    }
    return hr;
}

static HRESULT DumpShellMetadata(IMFMediaSource *pSource)
{
    HRESULT hr = S_OK;
    IMFGetService* pGetService = NULL;
    IPropertyStore *pProps = NULL;

    DebuggerPrint("--- Shell metadata ---\n");
    hr = pSource->QueryInterface(&pGetService);
    if (FAILED(hr))
    {
        DebuggerPrint("Media source doesn't support GetService.\n");
        goto Exit;
    }
    hr = pGetService->GetService(MF_PROPERTY_HANDLER_SERVICE, IID_PPV_ARGS(&pProps));
    if (FAILED(hr))
    {
        DebuggerPrint("Media source doesn't support property handler service.\n");
        goto Exit;
    }

    DWORD cProps;
    VT_HR_EXIT( pProps->GetCount(&cProps) );

    for (DWORD i = 0; i < cProps; i++)
    {
        PROPERTYKEY key;
        VT_HR_EXIT( pProps->GetAt(i, &key) );

        PROPVARIANT pv;
        VT_HR_EXIT( pProps->GetValue(key, &pv) );

        DisplayProperty(key, pv);
        PropVariantClear(&pv);
    }

Exit:
    SafeRelease(&pGetService);
    SafeRelease(&pProps);
    return hr;
}

static HRESULT DumpMetadata(const WCHAR* wszFileName)
{
    HRESULT hr = S_OK;
    IMFMediaSource* pSource = NULL;

    VT_HR_EXIT( CreateMediaSource(wszFileName, &pSource, NULL) );

    DebuggerPrint("\nMetadata for %S\n", wszFileName);
    DumpMediaMetadata(pSource);
    DumpShellMetadata(pSource);
    DebuggerPrint("\n");

Exit:
    SafeRelease(&pSource);
    return hr;
}

#endif // DEBUG_METADATA

#ifdef DEBUG_LOG_ATTRIBUTES

HRESULT LogAllMediaTypes(IMFSourceReader* pReader)
{
    HRESULT hr = S_OK;

    for (DWORD i = 0; ; i++)
    {
        for (DWORD j = 0; ; j++)
        {
            IMFMediaType *pType = NULL;

            // Get the j'th media type for the i'th stream.
            hr = pReader->GetNativeMediaType(i, j, &pType);
            if (FAILED(hr))
            {
                break;
            }

            // Examine the media type.
            DBGMSG(L"\n==> Stream %d, media type %d:\n", i, j);
            hr = LogAttributes(pType);

            pType->Release();
        }
        if (hr == MF_E_NO_MORE_TYPES)
        {
            hr = S_OK; // No more types for this stream. Continue.
        }
        else if (hr == MF_E_INVALIDSTREAMNUMBER)
        {
            hr = S_OK;
            break;   // No more streams. Quit the loop.
        }
        else if (FAILED(hr))
        {
             break;  // Some other error.
        }
    }

    return hr;
}

HRESULT LogAttributes(IMFAttributes *pAttributes)
{
    UINT32 count = 0;

    HRESULT hr = pAttributes->GetCount(&count);
    if (FAILED(hr))
    {
        return hr;
    }

    if (count == 0)
    {
        DBGMSG(L"\tNo attributes.\n");
    }

    for (UINT32 i = 0; i < count; i++)
    {
        hr = LogAttributeValueByIndex(pAttributes, i);
        if (FAILED(hr))
        {
            break;
        }
    }
    return hr;
}

HRESULT LogAttributeValueByIndex(IMFAttributes *pAttr, DWORD index)
{
    WCHAR *pGuidName = NULL;
    WCHAR *pGuidValName = NULL;

    GUID guid = { 0 };

    PROPVARIANT var;
    PropVariantInit(&var);

    HRESULT hr = pAttr->GetItemByIndex(index, &guid, &var);
    if (FAILED(hr))
    {
        goto done;
    }

    hr = GetGUIDName(guid, &pGuidName);
    if (FAILED(hr))
    {
        goto done;
    }

    DBGMSG(L"\t%s\t", pGuidName);

    hr = SpecialCaseAttributeValue(guid, var);
    if (FAILED(hr))
    {
        goto done;
    }
    if (hr == S_FALSE)
    {
        switch (var.vt)
        {
        case VT_UI4:
            DBGMSG(L"%d", var.ulVal);
            break;

        case VT_UI8:
            DBGMSG(L"%I64d", var.uhVal);
            break;

        case VT_R8:
            DBGMSG(L"%f", var.dblVal);
            break;

        case VT_CLSID:
            hr = GetGUIDName(*var.puuid, &pGuidValName);
            if (SUCCEEDED(hr))
            {
                DBGMSG(pGuidValName);
            }
            break;

        case VT_LPWSTR:
            DBGMSG(var.pwszVal);
            break;

        case VT_VECTOR | VT_UI1:
            DBGMSG(L"<<byte array>>");
            break;

        case VT_UNKNOWN:
            DBGMSG(L"IUnknown");
            break;

        default:
            DBGMSG(L"Unexpected attribute type (vt = %d)", var.vt);
            break;
        }
    }

done:
    DBGMSG(L"\n");
    CoTaskMemFree(pGuidName);
    CoTaskMemFree(pGuidValName);
    PropVariantClear(&var);
    return hr;
}

HRESULT GetGUIDName(const GUID& guid, WCHAR **ppwsz)
{
    HRESULT hr = S_OK;
    WCHAR *pName = NULL;

    LPCWSTR pcwsz = GetGUIDNameConst(guid);
    if (pcwsz)
    {
        size_t cchLength = 0;
    
        hr = StringCchLength(pcwsz, STRSAFE_MAX_CCH, &cchLength);
        if (FAILED(hr))
        {
            goto done;
        }
        
        pName = (WCHAR*)CoTaskMemAlloc((cchLength + 1) * sizeof(WCHAR));

        if (pName == NULL)
        {
            hr = E_OUTOFMEMORY;
            goto done;
        }

        hr = StringCchCopy(pName, cchLength + 1, pcwsz);
        if (FAILED(hr))
        {
            goto done;
        }
    }
    else
    {
        hr = StringFromCLSID(guid, &pName);
    }

done:
    if (FAILED(hr))
    {
        *ppwsz = NULL;
        CoTaskMemFree(pName);
    }
    else
    {
        *ppwsz = pName;
    }
    return hr;
}

void LogUINT32AsUINT64(const PROPVARIANT& var)
{
    UINT32 uHigh = 0, uLow = 0;
    Unpack2UINT32AsUINT64(var.uhVal.QuadPart, &uHigh, &uLow);
    DBGMSG(L"%d x %d", uHigh, uLow);
}

float OffsetToFloat(const MFOffset& offset)
{
    return offset.value + (static_cast<float>(offset.fract) / 65536.0f);
}

HRESULT LogVideoArea(const PROPVARIANT& var)
{
    if (var.caub.cElems < sizeof(MFVideoArea))
    {
        return MF_E_BUFFERTOOSMALL;
    }

    MFVideoArea *pArea = (MFVideoArea*)var.caub.pElems;

    DBGMSG(L"(%f, %f) (%d, %d)", 
           OffsetToFloat(pArea->OffsetX), OffsetToFloat(pArea->OffsetY), 
           pArea->Area.cx, pArea->Area.cy);
    return S_OK;
}

// Handle certain known special cases.
HRESULT SpecialCaseAttributeValue(GUID guid, const PROPVARIANT& var)
{
    if ((guid == MF_MT_FRAME_RATE) || (guid == MF_MT_FRAME_RATE_RANGE_MAX) ||
        (guid == MF_MT_FRAME_RATE_RANGE_MIN) || (guid == MF_MT_FRAME_SIZE) ||
        (guid == MF_MT_PIXEL_ASPECT_RATIO))
    {
        // Attributes that contain two packed 32-bit values.
        LogUINT32AsUINT64(var);
    }
    else if ((guid == MF_MT_GEOMETRIC_APERTURE) || 
             (guid == MF_MT_MINIMUM_DISPLAY_APERTURE) || 
             (guid == MF_MT_PAN_SCAN_APERTURE))
    {
        // Attributes that an MFVideoArea structure.
        return LogVideoArea(var);
    }
    else
    {
        return S_FALSE;
    }
    return S_OK;
}

#undef OutputDebugStringW

void DBGMSG(PCWSTR format, ...)
{
    va_list args;
    va_start(args, format);

    WCHAR msg[MAX_PATH];

    if (SUCCEEDED(StringCbVPrintf(msg, sizeof(msg), format, args)))
    {
        OutputDebugStringW(msg);
    }
}

#ifndef IF_EQUAL_RETURN
#define IF_EQUAL_RETURN(param, val) if(val == param) return L#val
#endif

LPCWSTR GetGUIDNameConst(const GUID& guid)
{
    IF_EQUAL_RETURN(guid, MF_MT_MAJOR_TYPE);
    IF_EQUAL_RETURN(guid, MF_MT_MAJOR_TYPE);
    IF_EQUAL_RETURN(guid, MF_MT_SUBTYPE);
    IF_EQUAL_RETURN(guid, MF_MT_ALL_SAMPLES_INDEPENDENT);
    IF_EQUAL_RETURN(guid, MF_MT_FIXED_SIZE_SAMPLES);
    IF_EQUAL_RETURN(guid, MF_MT_COMPRESSED);
    IF_EQUAL_RETURN(guid, MF_MT_SAMPLE_SIZE);
    IF_EQUAL_RETURN(guid, MF_MT_WRAPPED_TYPE);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_NUM_CHANNELS);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_SAMPLES_PER_SECOND);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_FLOAT_SAMPLES_PER_SECOND);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_AVG_BYTES_PER_SECOND);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_BLOCK_ALIGNMENT);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_BITS_PER_SAMPLE);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_VALID_BITS_PER_SAMPLE);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_SAMPLES_PER_BLOCK);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_CHANNEL_MASK);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_FOLDDOWN_MATRIX);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_WMADRC_PEAKREF);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_WMADRC_PEAKTARGET);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_WMADRC_AVGREF);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_WMADRC_AVGTARGET);
    IF_EQUAL_RETURN(guid, MF_MT_AUDIO_PREFER_WAVEFORMATEX);
    IF_EQUAL_RETURN(guid, MF_MT_AAC_PAYLOAD_TYPE);
    IF_EQUAL_RETURN(guid, MF_MT_AAC_AUDIO_PROFILE_LEVEL_INDICATION);
    IF_EQUAL_RETURN(guid, MF_MT_FRAME_SIZE);
    IF_EQUAL_RETURN(guid, MF_MT_FRAME_RATE);
    IF_EQUAL_RETURN(guid, MF_MT_FRAME_RATE_RANGE_MAX);
    IF_EQUAL_RETURN(guid, MF_MT_FRAME_RATE_RANGE_MIN);
    IF_EQUAL_RETURN(guid, MF_MT_PIXEL_ASPECT_RATIO);
    IF_EQUAL_RETURN(guid, MF_MT_DRM_FLAGS);
    IF_EQUAL_RETURN(guid, MF_MT_PAD_CONTROL_FLAGS);
    IF_EQUAL_RETURN(guid, MF_MT_SOURCE_CONTENT_HINT);
    IF_EQUAL_RETURN(guid, MF_MT_VIDEO_CHROMA_SITING);
    IF_EQUAL_RETURN(guid, MF_MT_INTERLACE_MODE);
    IF_EQUAL_RETURN(guid, MF_MT_TRANSFER_FUNCTION);
    IF_EQUAL_RETURN(guid, MF_MT_VIDEO_PRIMARIES);
    IF_EQUAL_RETURN(guid, MF_MT_CUSTOM_VIDEO_PRIMARIES);
    IF_EQUAL_RETURN(guid, MF_MT_YUV_MATRIX);
    IF_EQUAL_RETURN(guid, MF_MT_VIDEO_LIGHTING);
    IF_EQUAL_RETURN(guid, MF_MT_VIDEO_NOMINAL_RANGE);
    IF_EQUAL_RETURN(guid, MF_MT_GEOMETRIC_APERTURE);
    IF_EQUAL_RETURN(guid, MF_MT_MINIMUM_DISPLAY_APERTURE);
    IF_EQUAL_RETURN(guid, MF_MT_PAN_SCAN_APERTURE);
    IF_EQUAL_RETURN(guid, MF_MT_PAN_SCAN_ENABLED);
    IF_EQUAL_RETURN(guid, MF_MT_AVG_BITRATE);
    IF_EQUAL_RETURN(guid, MF_MT_AVG_BIT_ERROR_RATE);
    IF_EQUAL_RETURN(guid, MF_MT_MAX_KEYFRAME_SPACING);
    IF_EQUAL_RETURN(guid, MF_MT_DEFAULT_STRIDE);
    IF_EQUAL_RETURN(guid, MF_MT_PALETTE);
    IF_EQUAL_RETURN(guid, MF_MT_USER_DATA);
    IF_EQUAL_RETURN(guid, MF_MT_AM_FORMAT_TYPE);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG_START_TIME_CODE);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG2_PROFILE);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG2_LEVEL);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG2_FLAGS);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG_SEQUENCE_HEADER);
    IF_EQUAL_RETURN(guid, MF_MT_DV_AAUX_SRC_PACK_0);
    IF_EQUAL_RETURN(guid, MF_MT_DV_AAUX_CTRL_PACK_0);
    IF_EQUAL_RETURN(guid, MF_MT_DV_AAUX_SRC_PACK_1);
    IF_EQUAL_RETURN(guid, MF_MT_DV_AAUX_CTRL_PACK_1);
    IF_EQUAL_RETURN(guid, MF_MT_DV_VAUX_SRC_PACK);
    IF_EQUAL_RETURN(guid, MF_MT_DV_VAUX_CTRL_PACK);
    IF_EQUAL_RETURN(guid, MF_MT_ARBITRARY_HEADER);
    IF_EQUAL_RETURN(guid, MF_MT_ARBITRARY_FORMAT);
    IF_EQUAL_RETURN(guid, MF_MT_IMAGE_LOSS_TOLERANT); 
    IF_EQUAL_RETURN(guid, MF_MT_MPEG4_SAMPLE_DESCRIPTION);
    IF_EQUAL_RETURN(guid, MF_MT_MPEG4_CURRENT_SAMPLE_ENTRY);
    IF_EQUAL_RETURN(guid, MF_MT_ORIGINAL_4CC); 
    IF_EQUAL_RETURN(guid, MF_MT_ORIGINAL_WAVE_FORMAT_TAG);
    
    // Media types

    IF_EQUAL_RETURN(guid, MFMediaType_Audio);
    IF_EQUAL_RETURN(guid, MFMediaType_Video);
    IF_EQUAL_RETURN(guid, MFMediaType_Protected);
    IF_EQUAL_RETURN(guid, MFMediaType_SAMI);
    IF_EQUAL_RETURN(guid, MFMediaType_Script);
    IF_EQUAL_RETURN(guid, MFMediaType_Image);
    IF_EQUAL_RETURN(guid, MFMediaType_HTML);
    IF_EQUAL_RETURN(guid, MFMediaType_Binary);
    IF_EQUAL_RETURN(guid, MFMediaType_FileTransfer);

    IF_EQUAL_RETURN(guid, MFVideoFormat_RGB32); //    D3DFMT_X8R8G8B8 );
    IF_EQUAL_RETURN(guid, MFVideoFormat_ARGB32); //   D3DFMT_A8R8G8B8 );
    IF_EQUAL_RETURN(guid, MFVideoFormat_RGB24); //    D3DFMT_R8G8B8 );
    IF_EQUAL_RETURN(guid, MFVideoFormat_RGB555); //   D3DFMT_X1R5G5B5 );
    IF_EQUAL_RETURN(guid, MFVideoFormat_RGB565); //   D3DFMT_R5G6B5 );
    IF_EQUAL_RETURN(guid, MFVideoFormat_AI44); //     FCC('AI44') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_AYUV); //     FCC('AYUV') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_YUY2); //     FCC('YUY2') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_UYVY); //     FCC('UYVY') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_NV11); //     FCC('NV11') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_NV12); //     FCC('NV12') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_YV12); //     FCC('YV12') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_I420); //     FCC('I420') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_IYUV); //     FCC('IYUV') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_Y210); //     FCC('Y210') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_Y216); //     FCC('Y216') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_Y410); //     FCC('Y410') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_Y416); //     FCC('Y416') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_P210); //     FCC('P210') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_P216); //     FCC('P216') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_P010); //     FCC('P010') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_P016); //     FCC('P016') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_v210); //     FCC('v210') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_v410); //     FCC('v410') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MP43); //     FCC('MP43') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MP4S); //     FCC('MP4S') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_M4S2); //     FCC('M4S2') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MP4V); //     FCC('MP4V') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_WMV1); //     FCC('WMV1') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_WMV2); //     FCC('WMV2') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_WMV3); //     FCC('WMV3') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_WVC1); //     FCC('WVC1') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MSS1); //     FCC('MSS1') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MSS2); //     FCC('MSS2') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MPG1); //     FCC('MPG1') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_DVSL); //     FCC('dvsl') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_DVSD); //     FCC('dvsd') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_DV25); //     FCC('dv25') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_DV50); //     FCC('dv50') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_DVH1); //     FCC('dvh1') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_H264); //     FCC('H264') );
    IF_EQUAL_RETURN(guid, MFVideoFormat_MJPG);

    IF_EQUAL_RETURN(guid, MFAudioFormat_PCM); //              WAVE_FORMAT_PCM );
    IF_EQUAL_RETURN(guid, MFAudioFormat_Float); //            WAVE_FORMAT_IEEE_FLOAT );
    IF_EQUAL_RETURN(guid, MFAudioFormat_DTS); //              WAVE_FORMAT_DTS );
    IF_EQUAL_RETURN(guid, MFAudioFormat_Dolby_AC3_SPDIF); //  WAVE_FORMAT_DOLBY_AC3_SPDIF );
    IF_EQUAL_RETURN(guid, MFAudioFormat_DRM); //              WAVE_FORMAT_DRM );
    IF_EQUAL_RETURN(guid, MFAudioFormat_WMAudioV8); //        WAVE_FORMAT_WMAUDIO2 );
    IF_EQUAL_RETURN(guid, MFAudioFormat_WMAudioV9); //        WAVE_FORMAT_WMAUDIO3 );
    IF_EQUAL_RETURN(guid, MFAudioFormat_WMAudio_Lossless); // WAVE_FORMAT_WMAUDIO_LOSSLESS );
    IF_EQUAL_RETURN(guid, MFAudioFormat_WMASPDIF); //         WAVE_FORMAT_WMASPDIF );
    IF_EQUAL_RETURN(guid, MFAudioFormat_MSP1); //             WAVE_FORMAT_WMAVOICE9 );
    IF_EQUAL_RETURN(guid, MFAudioFormat_MP3); //              WAVE_FORMAT_MPEGLAYER3 );
    IF_EQUAL_RETURN(guid, MFAudioFormat_MPEG); //             WAVE_FORMAT_MPEG );
    IF_EQUAL_RETURN(guid, MFAudioFormat_AAC); //              WAVE_FORMAT_MPEG_HEAAC );
    IF_EQUAL_RETURN(guid, MFAudioFormat_ADTS); //             WAVE_FORMAT_MPEG_ADTS_AAC );

    return NULL;
}

#endif // DEBUG_LOG_ATTRIBUTES

CMFVideoDst::CMFVideoDst()
    : m_iRefCount(1)
    , m_bComStarted(false)
    , m_bMFInitialized(false)
    , m_pSinkWriter(NULL)
    , m_iWidth(0)
    , m_iHeight(0)
    , m_iFramesPerSecond(0)
    , m_streamIndex(0)
    , m_rtStart(0)
    , m_rtDuration(0)
    , m_iBitRate(0)
    , m_eVideoFormat(VFNone)
{
    HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
    m_bComStarted = SUCCEEDED(hr);
}

CMFVideoDst::~CMFVideoDst()
{
    Close();
    ShutdownMediaFoundation();
    if( m_bComStarted )
    {
        CoUninitialize();
    }
}

HRESULT CMFVideoDst::Clone(IVideoDst** ppClone)
{
    return vt::VtCreateVideoDst(ppClone);
}

HRESULT CMFVideoDst::OpenFile(__in_z const wchar_t* pwcFilename, int iValidPixelWidth, int iValidPixelHeight,
                              VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
                              double dAspectRatio, int iInterlaceMode)
{
    VT_HR_BEGIN();

    VT_HR_EXIT(InitializeMediaFoundation());

    VT_HR_EXIT(OpenByteStreamOrFile(NULL, pwcFilename, iValidPixelWidth, iValidPixelHeight,
        eVideoFormat, iFramesPerSecond, fBitsPerPixel, dAspectRatio, iInterlaceMode));

    VT_HR_EXIT_LABEL();

    return hr;
}

#if defined(VT_WINRT)
HRESULT CMFVideoDst::OpenFile(IStorageFile^ storageFile, int iValidPixelWidth, int iValidPixelHeight,
                              VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
                              double dAspectRatio, int iInterlaceMode)
{
    // The following code works well for WMV files.  Unfortunately, it creates MP4 files that
    // will play but don't contain a duration or bit rate.
#if 0
    IMFByteStream* pByteStream;

    VT_HR_BEGIN();

    VT_HR_EXIT(InitializeMediaFoundation());

    // Open the file for reading and writing.
    IRandomAccessStream^ stream = CSystem::Await(storageFile->OpenAsync(FileAccessMode::ReadWrite));

    // Create a byte stream from the random access stream.
    VT_HR_EXIT(MFCreateMFByteStreamOnStreamEx(reinterpret_cast<IUnknown *>(stream), &pByteStream));

    // Open the video from the byte stream.
    VT_HR_EXIT(OpenByteStreamOrFile(pByteStream, storageFile->Path->Data(), iValidPixelWidth, iValidPixelHeight,
        eVideoFormat, iFramesPerSecond, fBitsPerPixel, dAspectRatio, iInterlaceMode));

    VT_HR_EXIT_LABEL();

    if (FAILED(hr))
    {
        Close();
    }
    SafeRelease(&pByteStream);

    return hr;
#else
    VT_HR_BEGIN();

    VT_HR_EXIT(m_tempFileHelper.Initialize(storageFile, FileAccessMode::ReadWrite));

    IStorageFile^ tempStorageFile = m_tempFileHelper.GetTempStorageFile();

    VT_HR_EXIT(OpenFile(tempStorageFile->Path->Data(), iValidPixelWidth, iValidPixelHeight,
        eVideoFormat, iFramesPerSecond, fBitsPerPixel, dAspectRatio, iInterlaceMode));

    VT_HR_END();
#endif
}
#endif

HRESULT CMFVideoDst::WriteFrame(CImg &img)
{
    if (!m_pSinkWriter)
    {
        // Not opened yet.
        return E_UNEXPECTED;
    }
    if (m_eVideoFormat != CIMG)
    {
        return E_INVALIDARG;
    }
    HRESULT hr = S_OK;
    CRGBAImg imgWrappedBuffer;
    CRGBAImg imgOutput;
    CRGBAImg imgCropped;
    CRect rectCropRect;
    VT_HR_EXIT(VtConvertImage(imgOutput, img));
    if (m_iWidth > (UINT32) imgOutput.Width() || m_iHeight > (UINT32) imgOutput.Height())
    {
        return E_INVALIDARG;
    }        

    IMFSample *pSample = NULL;
    IMFMediaBuffer *pBuffer = NULL;

    const LONG cbWidth = 4 * m_iWidth;
    const DWORD cbBuffer = cbWidth * m_iHeight;

    BYTE *pData = NULL;

    // Create a new memory buffer.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateMemoryBuffer(cbBuffer, &pBuffer));

    // Lock the buffer and copy the video frame to the buffer.
    VT_HR_EXIT(pBuffer->Lock(&pData, NULL, NULL));
    VT_HR_EXIT(imgWrappedBuffer.Create(pData, m_iWidth, m_iHeight, cbWidth));
    rectCropRect = CRect(0, 0, m_iWidth, m_iHeight);
    VT_HR_EXIT(imgOutput.Share(imgCropped, &rectCropRect));
    VT_HR_EXIT(imgCropped.CopyTo(imgWrappedBuffer));
//    for(UINT32 iY=0; iY<m_iHeight; iY++)
//    {
//        VT_HR_EXIT(memcpy_s(pData, cbWidth, img.GetImg().Ptr(iY), cbWidth) == 0 ? S_OK : E_FAIL);
//        pData += cbWidth;
//    }

    if (pBuffer)
    {
        pBuffer->Unlock();
    }

    // Set the data length of the buffer.
    VT_HR_EXIT(pBuffer->SetCurrentLength(cbBuffer));

    // Create a media sample and add the buffer to the sample.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateSample(&pSample));
    VT_HR_EXIT(pSample->AddBuffer(pBuffer));

    // Set the time stamp and the duration.
    VT_HR_EXIT(pSample->SetSampleTime(m_rtStart));
    VT_HR_EXIT(pSample->SetSampleDuration(m_rtDuration));

    // Send the sample to the Sink Writer.
    VT_HR_EXIT(m_pSinkWriter->WriteSample(m_streamIndex, pSample));
    m_rtStart += m_rtDuration;

Exit:
    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
    return hr;
}

HRESULT CMFVideoDst::WriteFrame(CRGB32VideoImg &img)
{
    if (!m_pSinkWriter)
    {
        // Not opened yet.
        return E_UNEXPECTED;
    }
    // Work around for XVP color conversion problem
    // If try to send RGB32 to H264, the color converter is called
    // This seemingly uses the wrong color conversion leading to a shift in the output
    // In this case, 
    if (m_eVideoFormat == NV12 && m_videoOutputFormat == MFVideoFormat_H264)
    {
        CNV12VideoImg imgNV12;
        if (FAILED(VtConvertVideoImage(imgNV12, img)))
        {
            return E_FAIL;
        }
        else
        {
            return WriteFrame(imgNV12);
        }
    }
    if (m_eVideoFormat != RGB32)
    {
        return E_INVALIDARG;
    }
    CRect rectValidPixels = img.GetVideoImgInfo().rectValidPixels;
    if (rectValidPixels.left + m_iWidth > (UINT32) img.GetImg().Width() || rectValidPixels.top + m_iHeight > (UINT32) img.GetImg().Height())
    {
        return E_INVALIDARG;
    }        
    HRESULT hr = S_OK;
    CRGBAImg imgWrappedBuffer;
    CRGBAImg imgCropped;
    CRect rectCropRect;

    IMFSample *pSample = NULL;
    IMFMediaBuffer *pBuffer = NULL;

    const LONG cbWidth = 4 * m_iWidth;
    const DWORD cbBuffer = cbWidth * m_iHeight;

    BYTE *pData = NULL;

    // Create a new memory buffer.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateMemoryBuffer(cbBuffer, &pBuffer));

    // Lock the buffer and copy the video frame to the buffer.
    VT_HR_EXIT(pBuffer->Lock(&pData, NULL, NULL));
    VT_HR_EXIT(imgWrappedBuffer.Create(pData, m_iWidth, m_iHeight, cbWidth));
    rectCropRect = CRect(rectValidPixels.left, rectValidPixels.top, rectValidPixels.left+m_iWidth, rectValidPixels.top+m_iHeight);
    VT_HR_EXIT(img.GetImg().Share(imgCropped, &rectCropRect));
    VT_HR_EXIT(imgCropped.CopyTo(imgWrappedBuffer));
//    VT_HR_EXIT(img.GetImg().CopyTo(imgWrappedBuffer));
//    for(UINT32 iY=0; iY<m_iHeight; iY++)
//    {
//        VT_HR_EXIT(memcpy_s(pData, cbWidth, img.GetImg().Ptr(iY), cbWidth) == 0 ? S_OK : E_FAIL);
//        pData += cbWidth;
//    }

    if (pBuffer)
    {
        pBuffer->Unlock();
    }

    // Set the data length of the buffer.
    VT_HR_EXIT(pBuffer->SetCurrentLength(cbBuffer));

    // Create a media sample and add the buffer to the sample.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateSample(&pSample));
    VT_HR_EXIT(pSample->AddBuffer(pBuffer));

    // Set the time stamp and the duration.
    VT_HR_EXIT(pSample->SetSampleTime(m_rtStart));
    VT_HR_EXIT(pSample->SetSampleDuration(m_rtDuration));

    // Send the sample to the Sink Writer.
    VT_HR_EXIT(m_pSinkWriter->WriteSample(m_streamIndex, pSample));
    m_rtStart += m_rtDuration;

Exit:
    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
    return hr;
}

HRESULT CMFVideoDst::WriteFrame(CNV12VideoImg &img)
{
    if (!m_pSinkWriter)
    {
        // Not opened yet.
        return E_UNEXPECTED;
    }
    if (m_eVideoFormat != NV12)
    {
        return E_INVALIDARG;
    }
    CRect rectValidPixels = img.GetVideoImgInfo().rectValidPixels;
    if (rectValidPixels.left + m_iWidth > (UINT32) img.GetYImg().Width() || rectValidPixels.top + m_iHeight > (UINT32) img.GetYImg().Height())
    {
        return E_INVALIDARG;
    }        
    HRESULT hr = S_OK;
    CLumaByteImg imgWrappedBuffer;
    CUVByteImg imgWrappedUVBuffer;
    CLumaByteImg imgCroppedY;
    CUVByteImg imgCroppedUV;
    CRect rectCropRectY;
    CRect rectCropRectUV;

    IMFSample *pSample = NULL;
    IMFMediaBuffer *pBuffer = NULL;

    const LONG cbWidth = m_iWidth;
    LONG uvWidth = m_iWidth / 2;
    LONG uvHeight = m_iHeight / 2;
    const DWORD cbBuffer = cbWidth * m_iHeight + 2 * uvWidth * uvHeight;

    BYTE *pData = NULL;

    // Create a new memory buffer.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateMemoryBuffer(cbBuffer, &pBuffer));

    // Lock the buffer and copy the video frame to the buffer.
    VT_HR_EXIT(pBuffer->Lock(&pData, NULL, NULL));
    VT_HR_EXIT(imgWrappedBuffer.Create(pData, m_iWidth, m_iHeight, cbWidth));
    rectCropRectY = CRect(rectValidPixels.left, rectValidPixels.top, rectValidPixels.left+m_iWidth, rectValidPixels.top+m_iHeight);
    VT_HR_EXIT(img.GetYImg().Share(imgCroppedY, &rectCropRectY));
    VT_HR_EXIT(imgCroppedY.CopyTo(imgWrappedBuffer));
//    for(UINT32 iY=0; iY<m_iHeight; iY++)
//    {
//        VT_HR_EXIT(memcpy_s(pData, cbWidth, img.GetYImg().Ptr(iY), cbWidth) == 0 ? S_OK : E_FAIL);
//        pData += cbWidth;
//    }
    pData += m_iHeight*cbWidth;
    VT_HR_EXIT(imgWrappedUVBuffer.Create(pData, uvWidth, uvHeight, 2*uvWidth));
    rectCropRectUV = CRect(rectValidPixels.left/2, rectValidPixels.top/2, rectValidPixels.left/2+uvWidth, rectValidPixels.top/2+uvHeight);
    VT_HR_EXIT(img.GetUVImg().Share(imgCroppedUV, &rectCropRectUV));
    VT_HR_EXIT(imgCroppedUV.CopyTo(imgWrappedUVBuffer));
//    for(LONG iY=0; iY<uvHeight; iY++)
//    {
//        VT_HR_EXIT(memcpy_s(pData, 2*uvWidth, img.GetUVImg().Ptr(iY), 2*uvWidth) == 0 ? S_OK : E_FAIL);
//        pData += 2*uvWidth;
//    }

    if (pBuffer)
    {
        pBuffer->Unlock();
    }

    // Set the data length of the buffer.
    VT_HR_EXIT(pBuffer->SetCurrentLength(cbBuffer));

    // Create a media sample and add the buffer to the sample.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateSample(&pSample));
    VT_HR_EXIT(pSample->AddBuffer(pBuffer));

    // Set the time stamp and the duration.
    VT_HR_EXIT(pSample->SetSampleTime(m_rtStart));
    VT_HR_EXIT(pSample->SetSampleDuration(m_rtDuration));

    // Send the sample to the Sink Writer.
    VT_HR_EXIT(m_pSinkWriter->WriteSample(m_streamIndex, pSample));
    m_rtStart += m_rtDuration;

Exit:
    SafeRelease(&pSample);
    SafeRelease(&pBuffer);
    return hr;
}

HRESULT CMFVideoDst::Close()
{
    VT_HR_BEGIN();

    // Finish writing.
    if (m_pSinkWriter)
    {
        VT_HR_EXIT(m_pSinkWriter->Finalize());
    }

#if defined(VT_WINRT)
    VT_HR_EXIT(m_tempFileHelper.Finalize());
#endif

    VT_HR_EXIT_LABEL();

    SafeRelease(&m_pSinkWriter);
    m_iWidth = 0;
    m_iHeight = 0;
    m_iFramesPerSecond = 0;
    m_streamIndex = 0;
    m_rtStart = 0;
    m_rtDuration = 0;
    m_iBitRate = 0;
    m_eVideoFormat = VFNone;

    return hr;
}

ULONG CMFVideoDst::AddRef()
{
    return (ULONG)InterlockedIncrement(&m_iRefCount);
}

//-----------------------------------------------------------------------------
// Release
//-----------------------------------------------------------------------------

ULONG CMFVideoDst::Release()
{
    if(InterlockedDecrement(&m_iRefCount)==0)
    {
        delete this;
        return 0;
    }
    return (ULONG)m_iRefCount;
}

HRESULT CMFVideoDst::OpenByteStreamOrFile(__in_opt IMFByteStream* pByteStream, __in_z const wchar_t *pwcFilename,
                                          int iValidPixelWidth, int iValidPixelHeight,
                                          VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
                                          double dAspectRatio, int iInterlaceMode)
{
    IMFAttributes* pAttributes = NULL;
    IMFMediaType* pMediaTypeOut = NULL;
    IMFMediaType* pMediaTypeIn = NULL;

    // Check state and parameters.
    if (m_pSinkWriter)
    {
        // Already open.
        return E_UNEXPECTED;
    }
    if (pwcFilename == NULL || iValidPixelWidth < 2 || iValidPixelHeight < 2 || iFramesPerSecond < 0 || fBitsPerPixel < 0 || dAspectRatio <= 0)
    {
        return E_INVALIDARG;
    }
    if (eVideoFormat != CIMG && eVideoFormat != RGB32 && eVideoFormat != NV12)
    {
        return E_INVALIDARG;
    }

    VT_HR_BEGIN();

    m_iWidth = iValidPixelWidth - iValidPixelWidth % 2;
    m_iHeight = iValidPixelHeight - iValidPixelHeight % 2;
    m_iFramesPerSecond = iFramesPerSecond;
    m_iBitRate = int(fBitsPerPixel * iValidPixelWidth * iValidPixelHeight);
    m_eVideoFormat = eVideoFormat;
    m_rtStart = 0;
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFFrameRateToAverageTimePerFrame(m_iFramesPerSecond, 1, &m_rtDuration));

    size_t szLen = wcslen(pwcFilename);
    if (szLen > 5 && (pwcFilename[szLen-1] == 'v' || pwcFilename[szLen-1] == 'V')
                  && (pwcFilename[szLen-2] == 'm' || pwcFilename[szLen-2] == 'M')
                  && (pwcFilename[szLen-3] == 'w' || pwcFilename[szLen-3] == 'W')
                  && (pwcFilename[szLen-4] == '.'))
    {
        m_videoOutputFormat = MFVideoFormat_WVC1;
    }    
    else if (szLen > 5 && (pwcFilename[szLen-1] == 'f' || pwcFilename[szLen-1] == 'F')
                       && (pwcFilename[szLen-2] == 's' || pwcFilename[szLen-2] == 'S')
                       && (pwcFilename[szLen-3] == 'a' || pwcFilename[szLen-3] == 'A')
                       && (pwcFilename[szLen-4] == '.'))
    {
        m_videoOutputFormat = MFVideoFormat_WVC1;
    }
    else
    {
        m_videoOutputFormat = MFVideoFormat_H264;
    }

    // Work around for XVP color conversion problem
    // If try to send RGB32 to H264, the color converter is called
    // This seemingly uses the wrong color conversion leading to a shift in the output
    // In this case, 
    if (m_videoOutputFormat == MFVideoFormat_H264 && m_eVideoFormat == RGB32)
    {
        m_eVideoFormat = NV12;
    }

    // Create the sink writer.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateAttributes(&pAttributes, 10));
    VT_HR_EXIT(pAttributes->SetUINT32(MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS, 1));
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateSinkWriterFromURL(pwcFilename, pByteStream, pAttributes, &m_pSinkWriter));

    // Set the output media type.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateMediaType(&pMediaTypeOut));   
    VT_HR_EXIT(pMediaTypeOut->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video));     
    VT_HR_EXIT(pMediaTypeOut->SetGUID(MF_MT_SUBTYPE, m_videoOutputFormat));   
    VT_HR_EXIT(pMediaTypeOut->SetUINT32(MF_MT_AVG_BITRATE, m_iBitRate));   
    VT_HR_EXIT(pMediaTypeOut->SetUINT32(MF_MT_INTERLACE_MODE, iInterlaceMode));   
    VT_HR_EXIT(MFSetAttributeSize(pMediaTypeOut, MF_MT_FRAME_SIZE, m_iWidth, m_iHeight));   
    VT_HR_EXIT(MFSetAttributeRatio(pMediaTypeOut, MF_MT_FRAME_RATE, m_iFramesPerSecond, 1));   
    VT_HR_EXIT(MFSetAttributeRatio(pMediaTypeOut, MF_MT_PIXEL_ASPECT_RATIO, int(1000*dAspectRatio), 1000)); 
    VT_HR_EXIT(m_pSinkWriter->AddStream(pMediaTypeOut, &m_streamIndex));   

    // Set the input media type.
    VT_HR_EXIT(g_MediaFoundationAPI.pfnMFCreateMediaType(&pMediaTypeIn));   
    VT_HR_EXIT(pMediaTypeIn->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video));   
    switch(m_eVideoFormat)
    {
    case CIMG:
    case RGB32:
        VT_HR_EXIT(pMediaTypeIn->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_RGB32));     
        VT_HR_EXIT(pMediaTypeIn->SetUINT32(MF_MT_DEFAULT_STRIDE, 4*m_iWidth));   
        break;
    case NV12:
        VT_HR_EXIT(pMediaTypeIn->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_NV12));     
        VT_HR_EXIT(pMediaTypeIn->SetUINT32(MF_MT_DEFAULT_STRIDE, m_iWidth));   
        break;
    default:
        VT_HR_EXIT(E_INVALIDARG);
        break;
    };
    VT_HR_EXIT(pMediaTypeIn->SetUINT32(MF_MT_INTERLACE_MODE, iInterlaceMode));
    VT_HR_EXIT(MFSetAttributeSize(pMediaTypeIn, MF_MT_FRAME_SIZE, m_iWidth, m_iHeight));
    VT_HR_EXIT(MFSetAttributeRatio(pMediaTypeIn, MF_MT_FRAME_RATE, m_iFramesPerSecond, 1));   
    VT_HR_EXIT(MFSetAttributeRatio(pMediaTypeIn, MF_MT_PIXEL_ASPECT_RATIO, int(1000*dAspectRatio), 1000));   
    VT_HR_EXIT(m_pSinkWriter->SetInputMediaType(m_streamIndex, pMediaTypeIn, NULL));   

    // Tell the sink writer to start accepting data.
    VT_HR_EXIT(m_pSinkWriter->BeginWriting());

    VT_HR_EXIT_LABEL();

    if (FAILED(hr))
    {
        Close();
    }
    SafeRelease(&pMediaTypeOut);
    SafeRelease(&pMediaTypeIn);
    SafeRelease(&pAttributes);
    return hr;
}

HRESULT CMFVideoDst::InitializeMediaFoundation()
{
    HRESULT hr = S_OK;
    if (!m_bMFInitialized)
    {
        if (!VtIsMediaFoundationAvailable())
        {
            return E_NOINTERFACE;
        }
        hr = g_MediaFoundationAPI.pfnMFStartup(MF_VERSION, MFSTARTUP_FULL);
        m_bMFInitialized = true;
    }
    return hr;
}

HRESULT CMFVideoDst::ShutdownMediaFoundation()
{    
    HRESULT hr = S_OK;
    if (m_bMFInitialized)
    {    
        hr = g_MediaFoundationAPI.pfnMFShutdown();
        m_bMFInitialized = false;
    }
    return hr;
}

HRESULT vt::VtCreateVideoDst(IVideoDst **ppVideoDst) 
{
    HRESULT hr = S_OK;

    if (ppVideoDst == NULL)
    {
        return E_POINTER;
    }
    
    if (!vt::VtIsMediaFoundationAvailable())
    {
        return E_NOINTERFACE;
    }

    IVideoDst *videoDst = VT_NOTHROWNEW CMFVideoDst();

    *ppVideoDst = videoDst;
    return hr;
}
